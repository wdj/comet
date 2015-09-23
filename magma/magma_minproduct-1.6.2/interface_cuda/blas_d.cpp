/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
 
       @author Mark Gates
       @generated from blas_z.cpp normal z -> d, Fri Jan 30 19:00:20 2015
*/

#include <stdlib.h>
#include <stdio.h>

#include "magma_minproduct.h"
#include "error.h"

#define REAL

#ifdef HAVE_CUBLAS

// ----------------------------------------
// Convert MAGMA_minproduct constants to CUBLAS v1 constants, which are the same as lapack.
// These must be static to avoid conflict with CUBLAS v2 translators.
#define cublas_trans_const( magma_minproduct_const ) lapacke_trans_const( magma_minproduct_const )
#define cublas_side_const(  magma_minproduct_const ) lapacke_side_const(  magma_minproduct_const )
#define cublas_diag_const(  magma_minproduct_const ) lapacke_diag_const(  magma_minproduct_const )
#define cublas_uplo_const(  magma_minproduct_const ) lapacke_uplo_const(  magma_minproduct_const )


// ========================================
// Level 1 BLAS

// --------------------
/** Returns index of element of vector x having max. absolute value;
    i.e., max (infinity) norm.
    
    @param[in]
    n       Number of elements in vector x. n >= 0.
            
    @param[in]
    dx      DOUBLE_PRECISION array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).
            
    @param[in]
    incx    Stride between consecutive elements of dx. incx > 0.
            
    @ingroup magma_minproduct_dblas1
*/
extern "C" magma_minproduct_int_t
magma_minproduct_idamax(
    magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx )
{
    return cublasIdamax( n, dx, incx );
}

// --------------------
/** Returns index of element of vector x having min. absolute value.

    @param[in]
    n       Number of elements in vector x. n >= 0.

    @param[in]
    dx      DOUBLE_PRECISION array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx > 0.

    @ingroup magma_minproduct_dblas1
*/
extern "C" magma_minproduct_int_t
magma_minproduct_idamin(
    magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx )
{
    return cublasIdamin( n, dx, incx );
}

// --------------------
/** Returns the sum of absolute values of vector x; i.e., one norm.

    @param[in]
    n       Number of elements in vector x. n >= 0.

    @param[in]
    dx      DOUBLE_PRECISION array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx > 0.

    @ingroup magma_minproduct_dblas1
*/
extern "C" double
magma_minproduct_dasum(
    magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx )
{
    return cublasDasum( n, dx, incx );
}

// --------------------
/** Constant times a vector plus a vector; \f$ y = \alpha x + y \f$.

    @param[in]
    n       Number of elements in vectors x and y. n >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in]
    dx      DOUBLE_PRECISION array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in,out]
    dy      DOUBLE_PRECISION array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @ingroup magma_minproduct_dblas1
*/
extern "C" void
magma_minproduct_daxpy(
    magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDouble_ptr       dy, magma_minproduct_int_t incy )
{
    cublasDaxpy( n, alpha, dx, incx, dy, incy );
}

// --------------------
/** Copy vector x to vector y; \f$ y = x \f$.

    @param[in]
    n       Number of elements in vectors x and y. n >= 0.

    @param[in]
    dx      DOUBLE_PRECISION array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[out]
    dy      DOUBLE_PRECISION array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @ingroup magma_minproduct_dblas1
*/
extern "C" void
magma_minproduct_dcopy(
    magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDouble_ptr       dy, magma_minproduct_int_t incy )
{
    cublasDcopy( n, dx, incx, dy, incy );
}

// --------------------
/** Returns dot product of vectors x and y; \f$ x^H y \f$.

    @param[in]
    n       Number of elements in vector x and y. n >= 0.

    @param[in]
    dx      DOUBLE_PRECISION array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in]
    dy      DOUBLE_PRECISION array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @ingroup magma_minproduct_dblas1
*/
extern "C"
double magma_minproduct_ddot(
    magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDouble_const_ptr dy, magma_minproduct_int_t incy )
{
    return cublasDdot( n, dx, incx, dy, incy );
}

#ifdef COMPLEX
// --------------------
/** Returns dot product (unconjugated) of vectors x and y; \f$ x^T y \f$.

    @param[in]
    n       Number of elements in vector x and y. n >= 0.

    @param[in]
    dx      DOUBLE_PRECISION array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in]
    dy      DOUBLE_PRECISION array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @ingroup magma_minproduct_dblas1
*/
extern "C"
double magma_minproduct_ddot(
    magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDouble_const_ptr dy, magma_minproduct_int_t incy )
{
    return cublasDdot( n, dx, incx, dy, incy );
}
#endif // COMPLEX

// --------------------
/** Returns 2-norm of vector x. Avoids unnecesary over/underflow.

    @param[in]
    n       Number of elements in vector x and y. n >= 0.

    @param[in]
    dx      DOUBLE_PRECISION array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx > 0.

    @ingroup magma_minproduct_dblas1
*/
extern "C" double
magma_minproduct_dnrm2(
    magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx )
{
    return cublasDnrm2( n, dx, incx );
}

// --------------------
/** Apply Givens plane rotation, where cos (c) is real and sin (s) is real.

    @param[in]
    n       Number of elements in vector x and y. n >= 0.

    @param[in,out]
    dx      DOUBLE_PRECISION array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).
            On output, overwritten with c*x + s*y.

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in,out]
    dy      DOUBLE_PRECISION array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).
            On output, overwritten with -conj(s)*x + c*y.

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @param[in]
    c       double. cosine.

    @param[in]
    s       DOUBLE_PRECISION. sine. c and s define a rotation
            [ c         s ]  where c*c + s*conj(s) = 1.
            [ -conj(s)  c ]

    @ingroup magma_minproduct_dblas1
*/
extern "C" void
magma_minproduct_drot(
    magma_minproduct_int_t n,
    magma_minproductDouble_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDouble_ptr dy, magma_minproduct_int_t incy,
    double c, double s )
{
    cublasDrot( n, dx, incx, dy, incy, c, s );
}

#ifdef COMPLEX
// --------------------
/** Apply Givens plane rotation, where cos (c) and sin (s) are real.

    @param[in]
    n       Number of elements in vector x and y. n >= 0.

    @param[in,out]
    dx      DOUBLE_PRECISION array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).
            On output, overwritten with c*x + s*y.

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in,out]
    dy      DOUBLE_PRECISION array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).
            On output, overwritten with -conj(s)*x + c*y.

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @param[in]
    c       double. cosine.

    @param[in]
    s       double. sine. c and s define a rotation
            [  c  s ]  where c*c + s*s = 1.
            [ -s  c ]

    @ingroup magma_minproduct_dblas1
*/
extern "C" void
magma_minproduct_drot(
    magma_minproduct_int_t n,
    magma_minproductDouble_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDouble_ptr dy, magma_minproduct_int_t incy,
    double c, double s )
{
    cublasDrot( n, dx, incx, dy, incy, c, s );
}
#endif // COMPLEX

#ifdef REAL
// --------------------
/** Apply modified plane rotation.

    @ingroup magma_minproduct_dblas1
*/
extern "C" void
magma_minproduct_drotm(
    magma_minproduct_int_t n,
    double *dx, magma_minproduct_int_t incx,
    double *dy, magma_minproduct_int_t incy,
    const double *param )
{
    cublasDrotm( n, dx, incx, dy, incy, param );
}

// --------------------
/** Generate modified plane rotation.

    @ingroup magma_minproduct_dblas1
*/
extern "C" void
magma_minproduct_drotmg(
    double *d1, double       *d2,
    double *x1, const double *y1,
    double *param )
{
    cublasDrotmg( d1, d2, x1, y1, param );
}
#endif // REAL

// --------------------
/** Scales a vector by a constant; \f$ x = \alpha x \f$.

    @param[in]
    n       Number of elements in vector x. n >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in,out]
    dx      DOUBLE_PRECISION array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx > 0.

    @ingroup magma_minproduct_dblas1
*/
extern "C" void
magma_minproduct_dscal(
    magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_ptr dx, magma_minproduct_int_t incx )
{
    cublasDscal( n, alpha, dx, incx );
}

#ifdef COMPLEX
// --------------------
/** Scales a vector by a real constant; \f$ x = \alpha x \f$.

    @param[in]
    n       Number of elements in vector x. n >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$ (real)

    @param[in,out]
    dx      DOUBLE_PRECISION array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx > 0.

    @ingroup magma_minproduct_dblas1
*/
extern "C" void
magma_minproduct_dscal(
    magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_ptr dx, magma_minproduct_int_t incx )
{
    cublasDscal( n, alpha, dx, incx );
}
#endif // COMPLEX

// --------------------
/** Swap vector x and y; \f$ x <-> y \f$.

    @param[in]
    n       Number of elements in vector x and y. n >= 0.

    @param[in,out]
    dx      DOUBLE_PRECISION array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in,out]
    dy      DOUBLE_PRECISION array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @ingroup magma_minproduct_dblas1
*/
extern "C" void
magma_minproduct_dswap(
    magma_minproduct_int_t n,
    magma_minproductDouble_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDouble_ptr dy, magma_minproduct_int_t incy )
{
    cublasDswap( n, dx, incx, dy, incy );
}


// ========================================
// Level 2 BLAS

// --------------------
/** Perform matrix-vector product.
        \f$ y = \alpha A   x + \beta y \f$  (transA == Magma_minproductNoTrans), or \n
        \f$ y = \alpha A^T x + \beta y \f$  (transA == Magma_minproductTrans),   or \n
        \f$ y = \alpha A^H x + \beta y \f$  (transA == Magma_minproductConjTrans).

    @param[in]
    transA  Operation to perform on A.

    @param[in]
    m       Number of rows of A. m >= 0.

    @param[in]
    n       Number of columns of A. n >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in]
    dA      DOUBLE_PRECISION array of dimension (ldda,n), ldda >= max(1,m).
            The m-by-n matrix A, on GPU device.

    @param[in]
    ldda    Leading dimension of dA.

    @param[in]
    dx      DOUBLE_PRECISION array on GPU device.
            If transA == Magma_minproductNoTrans, the n element vector x of dimension (1 + (n-1)*incx); \n
            otherwise,                 the m element vector x of dimension (1 + (m-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in]
    beta    Scalar \f$ \beta \f$

    @param[in,out]
    dy      DOUBLE_PRECISION array on GPU device.
            If transA == Magma_minproductNoTrans, the m element vector y of dimension (1 + (m-1)*incy); \n
            otherwise,                 the n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @ingroup magma_minproduct_dblas2
*/
extern "C" void
magma_minproduct_dgemv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx,
    double beta,
    magma_minproductDouble_ptr       dy, magma_minproduct_int_t incy )
{
    cublasDgemv(
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
    dx      DOUBLE_PRECISION array on GPU device.
            The m element vector x of dimension (1 + (m-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in]
    dy      DOUBLE_PRECISION array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @param[in,out]
    dA      DOUBLE_PRECISION array on GPU device.
            The m-by-n matrix A of dimension (ldda,n), ldda >= max(1,m).

    @param[in]
    ldda    Leading dimension of dA.

    @ingroup magma_minproduct_dblas2
*/
extern "C" void
magma_minproduct_dger(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDouble_const_ptr dy, magma_minproduct_int_t incy,
    magma_minproductDouble_ptr       dA, magma_minproduct_int_t ldda )
{
    cublasDger(
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
    dx      DOUBLE_PRECISION array on GPU device.
            The m element vector x of dimension (1 + (m-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in]
    dy      DOUBLE_PRECISION array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @param[in,out]
    dA      DOUBLE_PRECISION array of dimension (ldda,n), ldda >= max(1,m).
            The m-by-n matrix A, on GPU device.

    @param[in]
    ldda    Leading dimension of dA.

    @ingroup magma_minproduct_dblas2
*/
extern "C" void
magma_minproduct_dger(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDouble_const_ptr dy, magma_minproduct_int_t incy,
    magma_minproductDouble_ptr       dA, magma_minproduct_int_t ldda )
{
    cublasDger(
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
    dA      DOUBLE_PRECISION array of dimension (ldda,n), ldda >= max(1,n).
            The n-by-n matrix A, on GPU device.

    @param[in]
    ldda    Leading dimension of dA.

    @param[in]
    dx      DOUBLE_PRECISION array on GPU device.
            The m element vector x of dimension (1 + (m-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in]
    beta    Scalar \f$ \beta \f$

    @param[in,out]
    dy      DOUBLE_PRECISION array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @ingroup magma_minproduct_dblas2
*/
extern "C" void
magma_minproduct_dsymv(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx,
    double beta,
    magma_minproductDouble_ptr       dy, magma_minproduct_int_t incy )
{
    cublasDsymv(
        cublas_uplo_const( uplo ),
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
    dx      DOUBLE_PRECISION array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in,out]
    dA      DOUBLE_PRECISION array of dimension (ldda,n), ldda >= max(1,n).
            The n-by-n matrix A, on GPU device.

    @param[in]
    ldda    Leading dimension of dA.

    @ingroup magma_minproduct_dblas2
*/
extern "C" void
magma_minproduct_dsyr(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDouble_ptr       dA, magma_minproduct_int_t ldda )
{
    cublasDsyr(
        cublas_uplo_const( uplo ),
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
    dx      DOUBLE_PRECISION array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @param[in]
    dy      DOUBLE_PRECISION array on GPU device.
            The n element vector y of dimension (1 + (n-1)*incy).

    @param[in]
    incy    Stride between consecutive elements of dy. incy != 0.

    @param[in,out]
    dA      DOUBLE_PRECISION array of dimension (ldda,n), ldda >= max(1,n).
            The n-by-n matrix A, on GPU device.

    @param[in]
    ldda    Leading dimension of dA.

    @ingroup magma_minproduct_dblas2
*/
extern "C" void
magma_minproduct_dsyr2(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dx, magma_minproduct_int_t incx,
    magma_minproductDouble_const_ptr dy, magma_minproduct_int_t incy,
    magma_minproductDouble_ptr       dA, magma_minproduct_int_t ldda )
{
    cublasDsyr2(
        cublas_uplo_const( uplo ),
        n,
        alpha, dx, incx,
               dy, incy,
               dA, ldda );
}

// --------------------
/** Perform triangular matrix-vector product.
        \f$ x = A   x \f$  (trans == Magma_minproductNoTrans), or \n
        \f$ x = A^T x \f$  (trans == Magma_minproductTrans),   or \n
        \f$ x = A^H x \f$  (trans == Magma_minproductConjTrans).

    @param[in]
    uplo    Whether the upper or lower triangle of A is referenced. 

    @param[in]
    trans   Operation to perform on A.

    @param[in]
    diag    Whether the diagonal of A is assumed to be unit or non-unit.

    @param[in]
    n       Number of rows and columns of A. n >= 0.

    @param[in]
    dA      DOUBLE_PRECISION array of dimension (ldda,n), ldda >= max(1,n).
            The n-by-n matrix A, on GPU device.

    @param[in]
    ldda    Leading dimension of dA.

    @param[in]
    dx      DOUBLE_PRECISION array on GPU device.
            The n element vector x of dimension (1 + (n-1)*incx).

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @ingroup magma_minproduct_dblas2
*/
extern "C" void
magma_minproduct_dtrmv(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_diag_t diag,
    magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr       dx, magma_minproduct_int_t incx )
{
    cublasDtrmv(
        cublas_uplo_const( uplo ),
        cublas_trans_const( trans ),
        cublas_diag_const( diag ),
        n,
        dA, ldda,
        dx, incx );
}

// --------------------
/** Solve triangular matrix-vector system (one right-hand side).
        \f$ A   x = b \f$  (trans == Magma_minproductNoTrans), or \n
        \f$ A^T x = b \f$  (trans == Magma_minproductTrans),   or \n
        \f$ A^H x = b \f$  (trans == Magma_minproductConjTrans).

    @param[in]
    uplo    Whether the upper or lower triangle of A is referenced. 

    @param[in]
    trans   Operation to perform on A.

    @param[in]
    diag    Whether the diagonal of A is assumed to be unit or non-unit.

    @param[in]
    n       Number of rows and columns of A. n >= 0.

    @param[in]
    dA      DOUBLE_PRECISION array of dimension (ldda,n), ldda >= max(1,n).
            The n-by-n matrix A, on GPU device.

    @param[in]
    ldda    Leading dimension of dA.

    @param[in,out]
    dx      DOUBLE_PRECISION array on GPU device.
            On entry, the n element RHS vector b of dimension (1 + (n-1)*incx).
            On exit, overwritten with the solution vector x.

    @param[in]
    incx    Stride between consecutive elements of dx. incx != 0.

    @ingroup magma_minproduct_dblas2
*/
extern "C" void
magma_minproduct_dtrsv(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_diag_t diag,
    magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr       dx, magma_minproduct_int_t incx )
{
    cublasDtrsv(
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
    dA      DOUBLE_PRECISION array on GPU device.
            If transA == Magma_minproductNoTrans, the m-by-k matrix A of dimension (ldda,k), ldda >= max(1,m); \n
            otherwise,                 the k-by-m matrix A of dimension (ldda,m), ldda >= max(1,k).

    @param[in]
    ldda    Leading dimension of dA.

    @param[in]
    dB      DOUBLE_PRECISION array on GPU device.
            If transB == Magma_minproductNoTrans, the k-by-n matrix B of dimension (lddb,n), lddb >= max(1,k); \n
            otherwise,                 the n-by-k matrix B of dimension (lddb,k), lddb >= max(1,n).

    @param[in]
    lddb    Leading dimension of dB.

    @param[in]
    beta    Scalar \f$ \beta \f$

    @param[in,out]
    dC      DOUBLE_PRECISION array on GPU device.
            The m-by-n matrix C of dimension (lddc,n), lddc >= max(1,m).

    @param[in]
    lddc    Leading dimension of dC.

    @ingroup magma_minproduct_dblas3
*/
extern "C" void
magma_minproduct_dgemm(
    magma_minproduct_trans_t transA, magma_minproduct_trans_t transB,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_const_ptr dB, magma_minproduct_int_t lddb,
    double beta,
    magma_minproductDouble_ptr       dC, magma_minproduct_int_t lddc )
{
    cublasDgemm(
        cublas_trans_const( transA ),
        cublas_trans_const( transB ),
        m, n, k,
        alpha, dA, ldda,
               dB, lddb,
        beta,  dC, lddc );
}

// --------------------
/** Perform symmetric matrix-matrix product.
        \f$ C = \alpha A B + \beta C \f$ (side == Magma_minproductLeft), or \n
        \f$ C = \alpha B A + \beta C \f$ (side == Magma_minproductRight),   \n
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
    dA      DOUBLE_PRECISION array on GPU device.
            If side == Magma_minproductLeft, the m-by-m symmetric matrix A of dimension (ldda,m), ldda >= max(1,m); \n
            otherwise,            the n-by-n symmetric matrix A of dimension (ldda,n), ldda >= max(1,n).

    @param[in]
    ldda    Leading dimension of dA.

    @param[in]
    dB      DOUBLE_PRECISION array on GPU device.
            The m-by-n matrix B of dimension (lddb,n), lddb >= max(1,m).

    @param[in]
    lddb    Leading dimension of dB.

    @param[in]
    beta    Scalar \f$ \beta \f$

    @param[in,out]
    dC      DOUBLE_PRECISION array on GPU device.
            The m-by-n matrix C of dimension (lddc,n), lddc >= max(1,m).

    @param[in]
    lddc    Leading dimension of dC.

    @ingroup magma_minproduct_dblas3
*/
extern "C" void
magma_minproduct_dsymm(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_const_ptr dB, magma_minproduct_int_t lddb,
    double beta,
    magma_minproductDouble_ptr       dC, magma_minproduct_int_t lddc )
{
    cublasDsymm(
        cublas_side_const( side ),
        cublas_uplo_const( uplo ),
        m, n,
        alpha, dA, ldda,
               dB, lddb,
        beta,  dC, lddc );
}

// --------------------
/** Perform symmetric rank-k update.
        \f$ C = \alpha A A^T + \beta C \f$ (trans == Magma_minproductNoTrans), or \n
        \f$ C = \alpha A^T A + \beta C \f$ (trans == Magma_minproductTrans),      \n
        where \f$ C \f$ is symmetric.

    @param[in]
    uplo    Whether the upper or lower triangle of C is referenced.

    @param[in]
    trans   Operation to perform on A.

    @param[in]
    n       Number of rows and columns of C. n >= 0.

    @param[in]
    k       Number of columns of A (for Magma_minproductNoTrans) or rows of A (for Magma_minproductTrans). k >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in]
    dA      DOUBLE_PRECISION array on GPU device.
            If trans == Magma_minproductNoTrans, the n-by-k matrix A of dimension (ldda,k), ldda >= max(1,n); \n
            otherwise,                the k-by-n matrix A of dimension (ldda,n), ldda >= max(1,k).

    @param[in]
    ldda    Leading dimension of dA.

    @param[in]
    beta    Scalar \f$ \beta \f$

    @param[in,out]
    dC      DOUBLE_PRECISION array on GPU device.
            The n-by-n symmetric matrix C of dimension (lddc,n), lddc >= max(1,n).

    @param[in]
    lddc    Leading dimension of dC.

    @ingroup magma_minproduct_dblas3
*/
extern "C" void
magma_minproduct_dsyrk(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    double beta,
    magma_minproductDouble_ptr       dC, magma_minproduct_int_t lddc )
{
    cublasDsyrk(
        cublas_uplo_const( uplo ),
        cublas_trans_const( trans ),
        n, k,
        alpha, dA, ldda,
        beta,  dC, lddc );
}

// --------------------
/** Perform symmetric rank-2k update.
        \f$ C = \alpha A B^T + \alpha B A^T \beta C \f$ (trans == Magma_minproductNoTrans), or \n
        \f$ C = \alpha A^T B + \alpha B^T A \beta C \f$ (trans == Magma_minproductTrans),      \n
        where \f$ C \f$ is symmetric.

    @param[in]
    uplo    Whether the upper or lower triangle of C is referenced.

    @param[in]
    trans   Operation to perform on A and B.

    @param[in]
    n       Number of rows and columns of C. n >= 0.

    @param[in]
    k       Number of columns of A and B (for Magma_minproductNoTrans) or rows of A and B (for Magma_minproductTrans). k >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in]
    dA      DOUBLE_PRECISION array on GPU device.
            If trans == Magma_minproductNoTrans, the n-by-k matrix A of dimension (ldda,k), ldda >= max(1,n); \n
            otherwise,                the k-by-n matrix A of dimension (ldda,n), ldda >= max(1,k).

    @param[in]
    ldda    Leading dimension of dA.

    @param[in]
    dB      DOUBLE_PRECISION array on GPU device.
            If trans == Magma_minproductNoTrans, the n-by-k matrix B of dimension (lddb,k), lddb >= max(1,n); \n
            otherwise,                the k-by-n matrix B of dimension (lddb,n), lddb >= max(1,k).

    @param[in]
    lddb    Leading dimension of dB.

    @param[in]
    beta    Scalar \f$ \beta \f$

    @param[in,out]
    dC      DOUBLE_PRECISION array on GPU device.
            The n-by-n symmetric matrix C of dimension (lddc,n), lddc >= max(1,n).

    @param[in]
    lddc    Leading dimension of dC.

    @ingroup magma_minproduct_dblas3
*/
extern "C" void
magma_minproduct_dsyr2k(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_const_ptr dB, magma_minproduct_int_t lddb,
    double beta,
    magma_minproductDouble_ptr       dC, magma_minproduct_int_t lddc )
{
    cublasDsyr2k(
        cublas_uplo_const( uplo ),
        cublas_trans_const( trans ),
        n, k,
        alpha, dA, ldda,
               dB, lddb,
        beta,  dC, lddc );
}

#ifdef COMPLEX
// --------------------
/** Perform symmetric matrix-matrix product.
        \f$ C = \alpha A B + \beta C \f$ (side == Magma_minproductLeft), or \n
        \f$ C = \alpha B A + \beta C \f$ (side == Magma_minproductRight),   \n
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
    dA      DOUBLE_PRECISION array on GPU device.
            If side == Magma_minproductLeft, the m-by-m symmetric matrix A of dimension (ldda,m), ldda >= max(1,m); \n
            otherwise,            the n-by-n symmetric matrix A of dimension (ldda,n), ldda >= max(1,n).

    @param[in]
    ldda    Leading dimension of dA.

    @param[in]
    dB      DOUBLE_PRECISION array on GPU device.
            The m-by-n matrix B of dimension (lddb,n), lddb >= max(1,m).

    @param[in]
    lddb    Leading dimension of dB.

    @param[in]
    beta    Scalar \f$ \beta \f$

    @param[in,out]
    dC      DOUBLE_PRECISION array on GPU device.
            The m-by-n matrix C of dimension (lddc,n), lddc >= max(1,m).

    @param[in]
    lddc    Leading dimension of dC.

    @ingroup magma_minproduct_dblas3
*/
extern "C" void
magma_minproduct_dsymm(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_const_ptr dB, magma_minproduct_int_t lddb,
    double beta,
    magma_minproductDouble_ptr       dC, magma_minproduct_int_t lddc )
{
    cublasDsymm(
        cublas_side_const( side ),
        cublas_uplo_const( uplo ),
        m, n,
        alpha, dA, ldda,
               dB, lddb,
        beta,  dC, lddc );
}

// --------------------
/** Perform symmetric rank-k update.
        \f$ C = \alpha A A^T + \beta C \f$ (trans == Magma_minproductNoTrans), or \n
        \f$ C = \alpha A^T A + \beta C \f$ (trans == Magma_minproductTrans),      \n
        where \f$ C \f$ is symmetric.

    @param[in]
    uplo    Whether the upper or lower triangle of C is referenced.

    @param[in]
    trans   Operation to perform on A.

    @param[in]
    n       Number of rows and columns of C. n >= 0.

    @param[in]
    k       Number of columns of A (for Magma_minproductNoTrans) or rows of A (for Magma_minproductTrans). k >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in]
    dA      DOUBLE_PRECISION array on GPU device.
            If trans == Magma_minproductNoTrans, the n-by-k matrix A of dimension (ldda,k), ldda >= max(1,n); \n
            otherwise,                the k-by-n matrix A of dimension (ldda,n), ldda >= max(1,k).

    @param[in]
    ldda    Leading dimension of dA.

    @param[in]
    beta    Scalar \f$ \beta \f$

    @param[in,out]
    dC      DOUBLE_PRECISION array on GPU device.
            The n-by-n symmetric matrix C of dimension (lddc,n), lddc >= max(1,n).

    @param[in]
    lddc    Leading dimension of dC.

    @ingroup magma_minproduct_dblas3
*/
extern "C" void
magma_minproduct_dsyrk(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    double beta,
    magma_minproductDouble_ptr       dC, magma_minproduct_int_t lddc )
{
    cublasDsyrk(
        cublas_uplo_const( uplo ),
        cublas_trans_const( trans ),
        n, k,
        alpha, dA, ldda,
        beta,  dC, lddc );
}

// --------------------
/** Perform symmetric rank-2k update.
        \f$ C = \alpha A B^T + \alpha B A^T \beta C \f$ (trans == Magma_minproductNoTrans), or \n
        \f$ C = \alpha A^T B + \alpha B^T A \beta C \f$ (trans == Magma_minproductTrans),      \n
        where \f$ C \f$ is symmetric.

    @param[in]
    uplo    Whether the upper or lower triangle of C is referenced.

    @param[in]
    trans   Operation to perform on A and B.

    @param[in]
    n       Number of rows and columns of C. n >= 0.

    @param[in]
    k       Number of columns of A and B (for Magma_minproductNoTrans) or rows of A and B (for Magma_minproductTrans). k >= 0.

    @param[in]
    alpha   Scalar \f$ \alpha \f$

    @param[in]
    dA      DOUBLE_PRECISION array on GPU device.
            If trans == Magma_minproductNoTrans, the n-by-k matrix A of dimension (ldda,k), ldda >= max(1,n); \n
            otherwise,                the k-by-n matrix A of dimension (ldda,n), ldda >= max(1,k).

    @param[in]
    ldda    Leading dimension of dA.

    @param[in]
    dB      DOUBLE_PRECISION array on GPU device.
            If trans == Magma_minproductNoTrans, the n-by-k matrix B of dimension (lddb,k), lddb >= max(1,n); \n
            otherwise,                the k-by-n matrix B of dimension (lddb,n), lddb >= max(1,k).

    @param[in]
    lddb    Leading dimension of dB.

    @param[in]
    beta    Scalar \f$ \beta \f$

    @param[in,out]
    dC      DOUBLE_PRECISION array on GPU device.
            The n-by-n symmetric matrix C of dimension (lddc,n), lddc >= max(1,n).

    @param[in]
    lddc    Leading dimension of dC.

    @ingroup magma_minproduct_dblas3
*/
extern "C" void
magma_minproduct_dsyr2k(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans,
    magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_const_ptr dB, magma_minproduct_int_t lddb,
    double beta,
    magma_minproductDouble_ptr       dC, magma_minproduct_int_t lddc )
{
    cublasDsyr2k(
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
        \f$ B = \alpha op(A) B \f$ (side == Magma_minproductLeft), or \n
        \f$ B = \alpha B op(A) \f$ (side == Magma_minproductRight),   \n
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
    dA      DOUBLE_PRECISION array on GPU device.
            If side == Magma_minproductLeft, the n-by-n triangular matrix A of dimension (ldda,n), ldda >= max(1,n); \n
            otherwise,            the m-by-m triangular matrix A of dimension (ldda,m), ldda >= max(1,m).

    @param[in]
    ldda    Leading dimension of dA.

    @param[in]
    dB      DOUBLE_PRECISION array on GPU device.
            The m-by-n matrix B of dimension (lddb,n), lddb >= max(1,m).

    @param[in]
    lddb    Leading dimension of dB.

    @ingroup magma_minproduct_dblas3
*/
extern "C" void
magma_minproduct_dtrmm(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_diag_t diag,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr       dB, magma_minproduct_int_t lddb )
{
    cublasDtrmm(
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
        \f$ op(A) X = \alpha B \f$ (side == Magma_minproductLeft), or \n
        \f$ X op(A) = \alpha B \f$ (side == Magma_minproductRight),   \n
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
    dA      DOUBLE_PRECISION array on GPU device.
            If side == Magma_minproductLeft, the m-by-m triangular matrix A of dimension (ldda,m), ldda >= max(1,m); \n
            otherwise,            the n-by-n triangular matrix A of dimension (ldda,n), ldda >= max(1,n).

    @param[in]
    ldda    Leading dimension of dA.

    @param[in,out]
    dB      DOUBLE_PRECISION array on GPU device.
            On entry, m-by-n matrix B of dimension (lddb,n), lddb >= max(1,m).
            On exit, overwritten with the solution matrix X.

    @param[in]
    lddb    Leading dimension of dB.

    @ingroup magma_minproduct_dblas3
*/
extern "C" void
magma_minproduct_dtrsm(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_diag_t diag,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr       dB, magma_minproduct_int_t lddb )
{
    cublasDtrsm(
        cublas_side_const( side ),
        cublas_uplo_const( uplo ),
        cublas_trans_const( trans ),
        cublas_diag_const( diag ),
        m, n,
        alpha, dA, ldda,
               dB, lddb );
}

#endif // HAVE_CUBLAS

#undef REAL
