/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
 
       @author Mathieu Faverge
       @author Mark Gates
*/

#ifndef MAGMA_minproduct_OPERATORS_H
#define MAGMA_minproduct_OPERATORS_H

// __host__ and __device__ are defined in CUDA headers.
#include "magma_minproduct.h"

/* names to match C++ std complex functions */
__host__ __device__ static inline double real(const magma_minproductDoubleComplex &x) { return MAGMA_minproduct_Z_REAL(x); }
__host__ __device__ static inline float  real(const magma_minproductFloatComplex  &x) { return MAGMA_minproduct_C_REAL(x); }
__host__ __device__ static inline double real(const double             &x) { return x; }
__host__ __device__ static inline float  real(const float              &x) { return x; }

__host__ __device__ static inline double imag(const magma_minproductDoubleComplex &x) { return MAGMA_minproduct_Z_IMAG(x); }
__host__ __device__ static inline float  imag(const magma_minproductFloatComplex  &x) { return MAGMA_minproduct_C_IMAG(x); }
__host__ __device__ static inline double imag(const double        & /*x*/) { return 0.; }
__host__ __device__ static inline float  imag(const float         & /*x*/) { return 0.; }

__host__ __device__ static inline magma_minproductDoubleComplex conj(const magma_minproductDoubleComplex &x) { return MAGMA_minproduct_Z_CNJG(x); }
__host__ __device__ static inline magma_minproductFloatComplex  conj(const magma_minproductFloatComplex  &x) { return MAGMA_minproduct_C_CNJG(x); }
__host__ __device__ static inline double             conj(const double             &x) { return x; }
__host__ __device__ static inline float              conj(const float              &x) { return x; }

__host__ __device__ static inline double fabs(const magma_minproductDoubleComplex &x) { return MAGMA_minproduct_Z_ABS(x); }
__host__ __device__ static inline float  fabs(const magma_minproductFloatComplex  &x) { return MAGMA_minproduct_C_ABS(x); }
//__host__ __device__ static inline float  fabs(const float              &x) { return MAGMA_minproduct_S_ABS(x); }  // conflicts with std::fabs

__host__ __device__ static inline double abs1(const magma_minproductDoubleComplex &x) { return MAGMA_minproduct_Z_ABS1(x); }
__host__ __device__ static inline float  abs1(const magma_minproductFloatComplex  &x) { return MAGMA_minproduct_C_ABS1(x); }
__host__ __device__ static inline double abs1(const double             &x) { return MAGMA_minproduct_D_ABS1(x); }
__host__ __device__ static inline float  abs1(const float              &x) { return MAGMA_minproduct_S_ABS1(x); }


/*************************************************************
 *              magma_minproductDoubleComplex
 */

// ---------- negate
__host__ __device__ static inline magma_minproductDoubleComplex
operator - (const magma_minproductDoubleComplex &a)
{
    return MAGMA_minproduct_Z_MAKE( -real(a),
                         -imag(a) );
}


// ---------- add
__host__ __device__ static inline magma_minproductDoubleComplex
operator + (const magma_minproductDoubleComplex a, const magma_minproductDoubleComplex b)
{
    return MAGMA_minproduct_Z_MAKE( real(a) + real(b),
                         imag(a) + imag(b) );
}

__host__ __device__ static inline magma_minproductDoubleComplex
operator + (const magma_minproductDoubleComplex a, const double s)
{
    return MAGMA_minproduct_Z_MAKE( real(a) + s,
                         imag(a) );
}

__host__ __device__ static inline magma_minproductDoubleComplex
operator + (const double s, const magma_minproductDoubleComplex b)
{
    return MAGMA_minproduct_Z_MAKE( s + real(b),
                             imag(b) );
}

__host__ __device__ static inline magma_minproductDoubleComplex&
operator += (magma_minproductDoubleComplex &a, const magma_minproductDoubleComplex b)
{
    a = MAGMA_minproduct_Z_MAKE( real(a) + real(b),
                      imag(a) + imag(b) );
    return a;
}

__host__ __device__ static inline magma_minproductDoubleComplex&
operator += (magma_minproductDoubleComplex &a, const double s)
{
    a = MAGMA_minproduct_Z_MAKE( real(a) + s,
                      imag(a) );
    return a;
}


// ---------- subtract
__host__ __device__ static inline magma_minproductDoubleComplex
operator - (const magma_minproductDoubleComplex a, const magma_minproductDoubleComplex b)
{
    return MAGMA_minproduct_Z_MAKE( real(a) - real(b),
                         imag(a) - imag(b) );
}

__host__ __device__ static inline magma_minproductDoubleComplex
operator - (const magma_minproductDoubleComplex a, const double s)
{
    return MAGMA_minproduct_Z_MAKE( real(a) - s,
                         imag(a) );
}

__host__ __device__ static inline magma_minproductDoubleComplex
operator - (const double s, const magma_minproductDoubleComplex b)
{
    return MAGMA_minproduct_Z_MAKE( s - real(b),
                           - imag(b) );
}

__host__ __device__ static inline magma_minproductDoubleComplex&
operator -= (magma_minproductDoubleComplex &a, const magma_minproductDoubleComplex b)
{
    a = MAGMA_minproduct_Z_MAKE( real(a) - real(b),
                      imag(a) - imag(b) );
    return a;
}

__host__ __device__ static inline magma_minproductDoubleComplex&
operator -= (magma_minproductDoubleComplex &a, const double s)
{
    a = MAGMA_minproduct_Z_MAKE( real(a) - s,
                      imag(a) );
    return a;
}


// ---------- multiply
__host__ __device__ static inline magma_minproductDoubleComplex
operator * (const magma_minproductDoubleComplex a, const magma_minproductDoubleComplex b)
{
    return MAGMA_minproduct_Z_MAKE( real(a)*real(b) - imag(a)*imag(b),
                         imag(a)*real(b) + real(a)*imag(b) );
}

__host__ __device__ static inline magma_minproductDoubleComplex
operator * (const magma_minproductDoubleComplex a, const double s)
{
    return MAGMA_minproduct_Z_MAKE( real(a)*s,
                         imag(a)*s );
}

__host__ __device__ static inline magma_minproductDoubleComplex
operator * (const double s, const magma_minproductDoubleComplex a)
{
    return MAGMA_minproduct_Z_MAKE( real(a)*s,
                         imag(a)*s );
}

__host__ __device__ static inline magma_minproductDoubleComplex&
operator *= (magma_minproductDoubleComplex &a, const magma_minproductDoubleComplex b)
{
    a = MAGMA_minproduct_Z_MAKE( real(a)*real(b) - imag(a)*imag(b),
                      imag(a)*real(b) + real(a)*imag(b) );
    return a;
}

__host__ __device__ static inline magma_minproductDoubleComplex&
operator *= (magma_minproductDoubleComplex &a, const double s)
{
    a = MAGMA_minproduct_Z_MAKE( real(a)*s,
                      imag(a)*s );
    return a;
}


// ---------- divide
/* From LAPACK DLADIV
 * Performs complex division in real arithmetic, avoiding unnecessary overflow.
 *
 *             a + i*b
 *  p + i*q = ---------
 *             c + i*d
 */
__host__ __device__ static inline magma_minproductDoubleComplex
operator / (const magma_minproductDoubleComplex x, const magma_minproductDoubleComplex y)
{
    double a = real(x);
    double b = imag(x);
    double c = real(y);
    double d = imag(y);
    double e, f, p, q;
    if ( fabs( d ) < fabs( c ) ) {
        e = d / c;
        f = c + d*e;
        p = ( a + b*e ) / f;
        q = ( b - a*e ) / f;
    }
    else {
        e = c / d;
        f = d + c*e;
        p = (  b + a*e ) / f;
        q = ( -a + b*e ) / f;
    }
    return MAGMA_minproduct_Z_MAKE( p, q );
}

__host__ __device__ static inline magma_minproductDoubleComplex
operator / (const magma_minproductDoubleComplex a, const double s)
{
    return MAGMA_minproduct_Z_MAKE( real(a)/s,
                         imag(a)/s );
}

__host__ __device__ static inline magma_minproductDoubleComplex
operator / (const double a, const magma_minproductDoubleComplex y)
{
    double c = real(y);
    double d = imag(y);
    double e, f, p, q;
    if ( fabs( d ) < fabs( c ) ) {
        e = d / c;
        f = c + d*e;
        p =  a   / f;
        q = -a*e / f;
    }
    else {
        e = c / d;
        f = d + c*e;
        p =  a*e / f;
        q = -a   / f;
    }
    return MAGMA_minproduct_Z_MAKE( p, q );
}

__host__ __device__ static inline magma_minproductDoubleComplex&
operator /= (magma_minproductDoubleComplex &a, const magma_minproductDoubleComplex b)
{
    a = a/b;
    return a;
}

__host__ __device__ static inline magma_minproductDoubleComplex&
operator /= (magma_minproductDoubleComplex &a, const double s)
{
    a = MAGMA_minproduct_Z_MAKE( real(a)/s,
                      imag(a)/s );
    return a;
}


// ---------- equality
__host__ __device__ static inline bool
operator == (const magma_minproductDoubleComplex a, const magma_minproductDoubleComplex b)
{
    return ( real(a) == real(b) &&
             imag(a) == imag(b) );
}

__host__ __device__ static inline bool
operator == (const magma_minproductDoubleComplex a, const double s)
{
    return ( real(a) == s &&
             imag(a) == 0. );
}

__host__ __device__ static inline bool
operator == (const double s, const magma_minproductDoubleComplex a)
{
    return ( real(a) == s &&
             imag(a) == 0. );
}


// ---------- not equality
__host__ __device__ static inline bool
operator != (const magma_minproductDoubleComplex a, const magma_minproductDoubleComplex b)
{
    return ! (a == b);
}

__host__ __device__ static inline bool
operator != (const magma_minproductDoubleComplex a, const double s)
{
    return ! (a == s);
}

__host__ __device__ static inline bool
operator != (const double s, const magma_minproductDoubleComplex a)
{
    return ! (a == s);
}


/*************************************************************
 *              magma_minproductFloatComplex
 */

// ---------- negate
__host__ __device__ static inline magma_minproductFloatComplex
operator - (const magma_minproductFloatComplex &a)
{
    return MAGMA_minproduct_C_MAKE( -real(a),
                         -imag(a) );
}


// ---------- add
__host__ __device__ static inline magma_minproductFloatComplex
operator + (const magma_minproductFloatComplex a, const magma_minproductFloatComplex b)
{
    return MAGMA_minproduct_C_MAKE( real(a) + real(b),
                         imag(a) + imag(b) );
}

__host__ __device__ static inline magma_minproductFloatComplex
operator + (const magma_minproductFloatComplex a, const float s)
{
    return MAGMA_minproduct_C_MAKE( real(a) + s,
                         imag(a) );
}

__host__ __device__ static inline magma_minproductFloatComplex
operator + (const float s, const magma_minproductFloatComplex b)
{
    return MAGMA_minproduct_C_MAKE( s + real(b),
                             imag(b) );
}

__host__ __device__ static inline magma_minproductFloatComplex&
operator += (magma_minproductFloatComplex &a, const magma_minproductFloatComplex b)
{
    a = MAGMA_minproduct_C_MAKE( real(a) + real(b),
                      imag(a) + imag(b) );
    return a;
}

__host__ __device__ static inline magma_minproductFloatComplex&
operator += (magma_minproductFloatComplex &a, const float s)
{
    a = MAGMA_minproduct_C_MAKE( real(a) + s,
                      imag(a) );
    return a;
}


// ---------- subtract
__host__ __device__ static inline magma_minproductFloatComplex
operator - (const magma_minproductFloatComplex a, const magma_minproductFloatComplex b)
{
    return MAGMA_minproduct_C_MAKE( real(a) - real(b),
                         imag(a) - imag(b) );
}

__host__ __device__ static inline magma_minproductFloatComplex
operator - (const magma_minproductFloatComplex a, const float s)
{
    return MAGMA_minproduct_C_MAKE( real(a) - s,
                         imag(a) );
}

__host__ __device__ static inline magma_minproductFloatComplex
operator - (const float s, const magma_minproductFloatComplex b)
{
    return MAGMA_minproduct_C_MAKE( s - real(b),
                           - imag(b) );
}

__host__ __device__ static inline magma_minproductFloatComplex&
operator -= (magma_minproductFloatComplex &a, const magma_minproductFloatComplex b)
{
    a = MAGMA_minproduct_C_MAKE( real(a) - real(b),
                      imag(a) - imag(b) );
    return a;
}

__host__ __device__ static inline magma_minproductFloatComplex&
operator -= (magma_minproductFloatComplex &a, const float s)
{
    a = MAGMA_minproduct_C_MAKE( real(a) - s,
                      imag(a) );
    return a;
}


// ---------- multiply
__host__ __device__ static inline magma_minproductFloatComplex
operator * (const magma_minproductFloatComplex a, const magma_minproductFloatComplex b)
{
    return MAGMA_minproduct_C_MAKE( real(a)*real(b) - imag(a)*imag(b),
                         imag(a)*real(b) + real(a)*imag(b) );
}

__host__ __device__ static inline magma_minproductFloatComplex
operator * (const magma_minproductFloatComplex a, const float s)
{
    return MAGMA_minproduct_C_MAKE( real(a)*s,
                         imag(a)*s );
}

__host__ __device__ static inline magma_minproductFloatComplex
operator * (const float s, const magma_minproductFloatComplex a)
{
    return MAGMA_minproduct_C_MAKE( real(a)*s,
                         imag(a)*s );
}

__host__ __device__ static inline magma_minproductFloatComplex&
operator *= (magma_minproductFloatComplex &a, const magma_minproductFloatComplex b)
{
    a = MAGMA_minproduct_C_MAKE( real(a)*real(b) - imag(a)*imag(b),
                      imag(a)*real(b) + real(a)*imag(b) );
    return a;
}

__host__ __device__ static inline magma_minproductFloatComplex&
operator *= (magma_minproductFloatComplex &a, const float s)
{
    a = MAGMA_minproduct_C_MAKE( real(a)*s,
                      imag(a)*s );
    return a;
}


// ---------- divide
/* From LAPACK DLADIV
 * Performs complex division in real arithmetic, avoiding unnecessary overflow.
 *
 *             a + i*b
 *  p + i*q = ---------
 *             c + i*d
 */
__host__ __device__ static inline magma_minproductFloatComplex
operator / (const magma_minproductFloatComplex x, const magma_minproductFloatComplex y)
{
    float a = real(x);
    float b = imag(x);
    float c = real(y);
    float d = imag(y);
    float e, f, p, q;
    if ( fabs( d ) < fabs( c ) ) {
        e = d / c;
        f = c + d*e;
        p = ( a + b*e ) / f;
        q = ( b - a*e ) / f;
    }
    else {
        e = c / d;
        f = d + c*e;
        p = (  b + a*e ) / f;
        q = ( -a + b*e ) / f;
    }
    return MAGMA_minproduct_C_MAKE( p, q );
}

__host__ __device__ static inline magma_minproductFloatComplex
operator / (const magma_minproductFloatComplex a, const float s)
{
    return MAGMA_minproduct_C_MAKE( real(a)/s,
                         imag(a)/s );
}

__host__ __device__ static inline magma_minproductFloatComplex
operator / (const float a, const magma_minproductFloatComplex y)
{
    float c = real(y);
    float d = imag(y);
    float e, f, p, q;
    if ( fabs( d ) < fabs( c ) ) {
        e = d / c;
        f = c + d*e;
        p =  a   / f;
        q = -a*e / f;
    }
    else {
        e = c / d;
        f = d + c*e;
        p =  a*e / f;
        q = -a   / f;
    }
    return MAGMA_minproduct_C_MAKE( p, q );
}

__host__ __device__ static inline magma_minproductFloatComplex&
operator /= (magma_minproductFloatComplex &a, const magma_minproductFloatComplex b)
{
    a = a/b;
    return a;
}

__host__ __device__ static inline magma_minproductFloatComplex&
operator /= (magma_minproductFloatComplex &a, const float s)
{
    a = MAGMA_minproduct_C_MAKE( real(a)/s,
                      imag(a)/s );
    return a;
}


// ---------- equality
__host__ __device__ static inline bool
operator == (const magma_minproductFloatComplex a, const magma_minproductFloatComplex b)
{
    return ( real(a) == real(b) &&
             imag(a) == imag(b) );
}

__host__ __device__ static inline bool
operator == (const magma_minproductFloatComplex a, const float s)
{
    return ( real(a) == s &&
             imag(a) == 0. );
}

__host__ __device__ static inline bool
operator == (const float s, const magma_minproductFloatComplex a)
{
    return ( real(a) == s &&
             imag(a) == 0. );
}


// ---------- not equality
__host__ __device__ static inline bool
operator != (const magma_minproductFloatComplex a, const magma_minproductFloatComplex b)
{
    return ! (a == b);
}

__host__ __device__ static inline bool
operator != (const magma_minproductFloatComplex a, const float s)
{
    return ! (a == s);
}

__host__ __device__ static inline bool
operator != (const float s, const magma_minproductFloatComplex a)
{
    return ! (a == s);
}

#endif /* MAGMA_minproduct_OPERATORS_H */
