/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
 
       @author Mathieu Faverge
       @author Mark Gates
*/

#ifndef MAGMA_tally4_OPERATORS_H
#define MAGMA_tally4_OPERATORS_H

// __host__ and __device__ are defined in CUDA headers.
#include "magma_tally4.h"

/* names to match C++ std complex functions */
__host__ __device__ static inline double real(const magma_tally4DoubleComplex &x) { return MAGMA_tally4_Z_REAL(x); }
__host__ __device__ static inline float  real(const magma_tally4FloatComplex  &x) { return MAGMA_tally4_C_REAL(x); }
__host__ __device__ static inline double real(const double             &x) { return x; }
__host__ __device__ static inline float  real(const float              &x) { return x; }

__host__ __device__ static inline double imag(const magma_tally4DoubleComplex &x) { return MAGMA_tally4_Z_IMAG(x); }
__host__ __device__ static inline float  imag(const magma_tally4FloatComplex  &x) { return MAGMA_tally4_C_IMAG(x); }
__host__ __device__ static inline double imag(const double        & /*x*/) { return 0.; }
__host__ __device__ static inline float  imag(const float         & /*x*/) { return 0.; }

__host__ __device__ static inline magma_tally4DoubleComplex conj(const magma_tally4DoubleComplex &x) { return MAGMA_tally4_Z_CNJG(x); }
__host__ __device__ static inline magma_tally4FloatComplex  conj(const magma_tally4FloatComplex  &x) { return MAGMA_tally4_C_CNJG(x); }
__host__ __device__ static inline double             conj(const double             &x) { return x; }
__host__ __device__ static inline float              conj(const float              &x) { return x; }

__host__ __device__ static inline double fabs(const magma_tally4DoubleComplex &x) { return MAGMA_tally4_Z_ABS(x); }
__host__ __device__ static inline float  fabs(const magma_tally4FloatComplex  &x) { return MAGMA_tally4_C_ABS(x); }
//__host__ __device__ static inline float  fabs(const float              &x) { return MAGMA_tally4_S_ABS(x); }  // conflicts with std::fabs

__host__ __device__ static inline double abs1(const magma_tally4DoubleComplex &x) { return MAGMA_tally4_Z_ABS1(x); }
__host__ __device__ static inline float  abs1(const magma_tally4FloatComplex  &x) { return MAGMA_tally4_C_ABS1(x); }
__host__ __device__ static inline double abs1(const double             &x) { return MAGMA_tally4_D_ABS1(x); }
__host__ __device__ static inline float  abs1(const float              &x) { return MAGMA_tally4_S_ABS1(x); }


/*************************************************************
 *              magma_tally4DoubleComplex
 */

// ---------- negate
__host__ __device__ static inline magma_tally4DoubleComplex
operator - (const magma_tally4DoubleComplex &a)
{
    return MAGMA_tally4_Z_MAKE( -real(a),
                         -imag(a) );
}


// ---------- add
__host__ __device__ static inline magma_tally4DoubleComplex
operator + (const magma_tally4DoubleComplex a, const magma_tally4DoubleComplex b)
{
    return MAGMA_tally4_Z_MAKE( real(a) + real(b),
                         imag(a) + imag(b) );
}

__host__ __device__ static inline magma_tally4DoubleComplex
operator + (const magma_tally4DoubleComplex a, const double s)
{
    return MAGMA_tally4_Z_MAKE( real(a) + s,
                         imag(a) );
}

__host__ __device__ static inline magma_tally4DoubleComplex
operator + (const double s, const magma_tally4DoubleComplex b)
{
    return MAGMA_tally4_Z_MAKE( s + real(b),
                             imag(b) );
}

__host__ __device__ static inline magma_tally4DoubleComplex&
operator += (magma_tally4DoubleComplex &a, const magma_tally4DoubleComplex b)
{
    a = MAGMA_tally4_Z_MAKE( real(a) + real(b),
                      imag(a) + imag(b) );
    return a;
}

__host__ __device__ static inline magma_tally4DoubleComplex&
operator += (magma_tally4DoubleComplex &a, const double s)
{
    a = MAGMA_tally4_Z_MAKE( real(a) + s,
                      imag(a) );
    return a;
}


// ---------- subtract
__host__ __device__ static inline magma_tally4DoubleComplex
operator - (const magma_tally4DoubleComplex a, const magma_tally4DoubleComplex b)
{
    return MAGMA_tally4_Z_MAKE( real(a) - real(b),
                         imag(a) - imag(b) );
}

__host__ __device__ static inline magma_tally4DoubleComplex
operator - (const magma_tally4DoubleComplex a, const double s)
{
    return MAGMA_tally4_Z_MAKE( real(a) - s,
                         imag(a) );
}

__host__ __device__ static inline magma_tally4DoubleComplex
operator - (const double s, const magma_tally4DoubleComplex b)
{
    return MAGMA_tally4_Z_MAKE( s - real(b),
                           - imag(b) );
}

__host__ __device__ static inline magma_tally4DoubleComplex&
operator -= (magma_tally4DoubleComplex &a, const magma_tally4DoubleComplex b)
{
    a = MAGMA_tally4_Z_MAKE( real(a) - real(b),
                      imag(a) - imag(b) );
    return a;
}

__host__ __device__ static inline magma_tally4DoubleComplex&
operator -= (magma_tally4DoubleComplex &a, const double s)
{
    a = MAGMA_tally4_Z_MAKE( real(a) - s,
                      imag(a) );
    return a;
}


// ---------- multiply
__host__ __device__ static inline magma_tally4DoubleComplex
operator * (const magma_tally4DoubleComplex a, const magma_tally4DoubleComplex b)
{
    return MAGMA_tally4_Z_MAKE( real(a)*real(b) - imag(a)*imag(b),
                         imag(a)*real(b) + real(a)*imag(b) );
}

__host__ __device__ static inline magma_tally4DoubleComplex
operator * (const magma_tally4DoubleComplex a, const double s)
{
    return MAGMA_tally4_Z_MAKE( real(a)*s,
                         imag(a)*s );
}

__host__ __device__ static inline magma_tally4DoubleComplex
operator * (const double s, const magma_tally4DoubleComplex a)
{
    return MAGMA_tally4_Z_MAKE( real(a)*s,
                         imag(a)*s );
}

__host__ __device__ static inline magma_tally4DoubleComplex&
operator *= (magma_tally4DoubleComplex &a, const magma_tally4DoubleComplex b)
{
    a = MAGMA_tally4_Z_MAKE( real(a)*real(b) - imag(a)*imag(b),
                      imag(a)*real(b) + real(a)*imag(b) );
    return a;
}

__host__ __device__ static inline magma_tally4DoubleComplex&
operator *= (magma_tally4DoubleComplex &a, const double s)
{
    a = MAGMA_tally4_Z_MAKE( real(a)*s,
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
__host__ __device__ static inline magma_tally4DoubleComplex
operator / (const magma_tally4DoubleComplex x, const magma_tally4DoubleComplex y)
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
    return MAGMA_tally4_Z_MAKE( p, q );
}

__host__ __device__ static inline magma_tally4DoubleComplex
operator / (const magma_tally4DoubleComplex a, const double s)
{
    return MAGMA_tally4_Z_MAKE( real(a)/s,
                         imag(a)/s );
}

__host__ __device__ static inline magma_tally4DoubleComplex
operator / (const double a, const magma_tally4DoubleComplex y)
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
    return MAGMA_tally4_Z_MAKE( p, q );
}

__host__ __device__ static inline magma_tally4DoubleComplex&
operator /= (magma_tally4DoubleComplex &a, const magma_tally4DoubleComplex b)
{
    a = a/b;
    return a;
}

__host__ __device__ static inline magma_tally4DoubleComplex&
operator /= (magma_tally4DoubleComplex &a, const double s)
{
    a = MAGMA_tally4_Z_MAKE( real(a)/s,
                      imag(a)/s );
    return a;
}


// ---------- equality
__host__ __device__ static inline bool
operator == (const magma_tally4DoubleComplex a, const magma_tally4DoubleComplex b)
{
    return ( real(a) == real(b) &&
             imag(a) == imag(b) );
}

__host__ __device__ static inline bool
operator == (const magma_tally4DoubleComplex a, const double s)
{
    return ( real(a) == s &&
             imag(a) == 0. );
}

__host__ __device__ static inline bool
operator == (const double s, const magma_tally4DoubleComplex a)
{
    return ( real(a) == s &&
             imag(a) == 0. );
}


// ---------- not equality
__host__ __device__ static inline bool
operator != (const magma_tally4DoubleComplex a, const magma_tally4DoubleComplex b)
{
    return ! (a == b);
}

__host__ __device__ static inline bool
operator != (const magma_tally4DoubleComplex a, const double s)
{
    return ! (a == s);
}

__host__ __device__ static inline bool
operator != (const double s, const magma_tally4DoubleComplex a)
{
    return ! (a == s);
}


/*************************************************************
 *              magma_tally4FloatComplex
 */

// ---------- negate
__host__ __device__ static inline magma_tally4FloatComplex
operator - (const magma_tally4FloatComplex &a)
{
    return MAGMA_tally4_C_MAKE( -real(a),
                         -imag(a) );
}


// ---------- add
__host__ __device__ static inline magma_tally4FloatComplex
operator + (const magma_tally4FloatComplex a, const magma_tally4FloatComplex b)
{
    return MAGMA_tally4_C_MAKE( real(a) + real(b),
                         imag(a) + imag(b) );
}

__host__ __device__ static inline magma_tally4FloatComplex
operator + (const magma_tally4FloatComplex a, const float s)
{
    return MAGMA_tally4_C_MAKE( real(a) + s,
                         imag(a) );
}

__host__ __device__ static inline magma_tally4FloatComplex
operator + (const float s, const magma_tally4FloatComplex b)
{
    return MAGMA_tally4_C_MAKE( s + real(b),
                             imag(b) );
}

__host__ __device__ static inline magma_tally4FloatComplex&
operator += (magma_tally4FloatComplex &a, const magma_tally4FloatComplex b)
{
    a = MAGMA_tally4_C_MAKE( real(a) + real(b),
                      imag(a) + imag(b) );
    return a;
}

__host__ __device__ static inline magma_tally4FloatComplex&
operator += (magma_tally4FloatComplex &a, const float s)
{
    a = MAGMA_tally4_C_MAKE( real(a) + s,
                      imag(a) );
    return a;
}


// ---------- subtract
__host__ __device__ static inline magma_tally4FloatComplex
operator - (const magma_tally4FloatComplex a, const magma_tally4FloatComplex b)
{
    return MAGMA_tally4_C_MAKE( real(a) - real(b),
                         imag(a) - imag(b) );
}

__host__ __device__ static inline magma_tally4FloatComplex
operator - (const magma_tally4FloatComplex a, const float s)
{
    return MAGMA_tally4_C_MAKE( real(a) - s,
                         imag(a) );
}

__host__ __device__ static inline magma_tally4FloatComplex
operator - (const float s, const magma_tally4FloatComplex b)
{
    return MAGMA_tally4_C_MAKE( s - real(b),
                           - imag(b) );
}

__host__ __device__ static inline magma_tally4FloatComplex&
operator -= (magma_tally4FloatComplex &a, const magma_tally4FloatComplex b)
{
    a = MAGMA_tally4_C_MAKE( real(a) - real(b),
                      imag(a) - imag(b) );
    return a;
}

__host__ __device__ static inline magma_tally4FloatComplex&
operator -= (magma_tally4FloatComplex &a, const float s)
{
    a = MAGMA_tally4_C_MAKE( real(a) - s,
                      imag(a) );
    return a;
}


// ---------- multiply
__host__ __device__ static inline magma_tally4FloatComplex
operator * (const magma_tally4FloatComplex a, const magma_tally4FloatComplex b)
{
    return MAGMA_tally4_C_MAKE( real(a)*real(b) - imag(a)*imag(b),
                         imag(a)*real(b) + real(a)*imag(b) );
}

__host__ __device__ static inline magma_tally4FloatComplex
operator * (const magma_tally4FloatComplex a, const float s)
{
    return MAGMA_tally4_C_MAKE( real(a)*s,
                         imag(a)*s );
}

__host__ __device__ static inline magma_tally4FloatComplex
operator * (const float s, const magma_tally4FloatComplex a)
{
    return MAGMA_tally4_C_MAKE( real(a)*s,
                         imag(a)*s );
}

__host__ __device__ static inline magma_tally4FloatComplex&
operator *= (magma_tally4FloatComplex &a, const magma_tally4FloatComplex b)
{
    a = MAGMA_tally4_C_MAKE( real(a)*real(b) - imag(a)*imag(b),
                      imag(a)*real(b) + real(a)*imag(b) );
    return a;
}

__host__ __device__ static inline magma_tally4FloatComplex&
operator *= (magma_tally4FloatComplex &a, const float s)
{
    a = MAGMA_tally4_C_MAKE( real(a)*s,
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
__host__ __device__ static inline magma_tally4FloatComplex
operator / (const magma_tally4FloatComplex x, const magma_tally4FloatComplex y)
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
    return MAGMA_tally4_C_MAKE( p, q );
}

__host__ __device__ static inline magma_tally4FloatComplex
operator / (const magma_tally4FloatComplex a, const float s)
{
    return MAGMA_tally4_C_MAKE( real(a)/s,
                         imag(a)/s );
}

__host__ __device__ static inline magma_tally4FloatComplex
operator / (const float a, const magma_tally4FloatComplex y)
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
    return MAGMA_tally4_C_MAKE( p, q );
}

__host__ __device__ static inline magma_tally4FloatComplex&
operator /= (magma_tally4FloatComplex &a, const magma_tally4FloatComplex b)
{
    a = a/b;
    return a;
}

__host__ __device__ static inline magma_tally4FloatComplex&
operator /= (magma_tally4FloatComplex &a, const float s)
{
    a = MAGMA_tally4_C_MAKE( real(a)/s,
                      imag(a)/s );
    return a;
}


// ---------- equality
__host__ __device__ static inline bool
operator == (const magma_tally4FloatComplex a, const magma_tally4FloatComplex b)
{
    return ( real(a) == real(b) &&
             imag(a) == imag(b) );
}

__host__ __device__ static inline bool
operator == (const magma_tally4FloatComplex a, const float s)
{
    return ( real(a) == s &&
             imag(a) == 0. );
}

__host__ __device__ static inline bool
operator == (const float s, const magma_tally4FloatComplex a)
{
    return ( real(a) == s &&
             imag(a) == 0. );
}


// ---------- not equality
__host__ __device__ static inline bool
operator != (const magma_tally4FloatComplex a, const magma_tally4FloatComplex b)
{
    return ! (a == b);
}

__host__ __device__ static inline bool
operator != (const magma_tally4FloatComplex a, const float s)
{
    return ! (a == s);
}

__host__ __device__ static inline bool
operator != (const float s, const magma_tally4FloatComplex a)
{
    return ! (a == s);
}

#endif /* MAGMA_tally4_OPERATORS_H */
