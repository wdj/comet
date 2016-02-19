/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015
*/

#ifndef MAGMA_tally3_TYPES_H
#define MAGMA_tally3_TYPES_H

#include <stdint.h>
#include <assert.h>


// each implementation of MAGMA_tally3 defines HAVE_* appropriately.
#if ! defined(HAVE_CUBLAS) && ! defined(HAVE_clAmdBlas) && ! defined(HAVE_MIC)
#define HAVE_CUBLAS
#endif


// ========================================
// C99 standard defines __func__. Some older compilers use __FUNCTION__.
// Note __func__ in C99 is not a macro, so ifndef __func__ doesn't work.
#if __STDC_VERSION__ < 199901L
  #ifndef __func__
    #if __GNUC__ >= 2 || _MSC_VER >= 1300
      #define __func__ __FUNCTION__
    #else
      #define __func__ "<unknown>"
    #endif
  #endif
#endif


// ========================================
// To use int64_t, link with mkl_intel_ilp64 or similar (instead of mkl_intel_lp64).
// Similar to magma_tally3_int_t we declare magma_tally3_index_t used for row/column indices in sparse
#if defined(MAGMA_tally3_ILP64) || defined(MKL_ILP64)
//typedef int64_t magma_tally3_int_t;
typedef long long int magma_tally3_int_t;  // MKL uses long long int, not int64_t
#else
typedef int magma_tally3_int_t;
#endif

typedef int magma_tally3_index_t;

// Define new type that the precision generator will not change (matches PLASMA)
typedef double real_Double_t;


// ========================================
// define types specific to implementation (CUDA, OpenCL, MIC)
// define macros to deal with complex numbers
#if defined(HAVE_CUBLAS)
    #ifndef CUBLAS_V2_H_
    #include <cublas.h>
    #endif
    
    typedef cudaStream_t   magma_tally3_queue_t;
    typedef cudaEvent_t    magma_tally3_event_t;
    typedef int            magma_tally3_device_t;
    
    typedef cuDoubleComplex magma_tally3DoubleComplex;
    typedef cuFloatComplex  magma_tally3FloatComplex;
    
    #define MAGMA_tally3_Z_MAKE(r,i)     make_cuDoubleComplex(r, i)
    #define MAGMA_tally3_Z_REAL(a)       (a).x
    #define MAGMA_tally3_Z_IMAG(a)       (a).y
    #define MAGMA_tally3_Z_ADD(a, b)     cuCadd(a, b)
    #define MAGMA_tally3_Z_SUB(a, b)     cuCsub(a, b)
    #define MAGMA_tally3_Z_MUL(a, b)     cuCmul(a, b)
    #define MAGMA_tally3_Z_DIV(a, b)     cuCdiv(a, b)
    #define MAGMA_tally3_Z_ABS(a)        cuCabs(a)
    #define MAGMA_tally3_Z_ABS1(a)       (fabs((a).x) + fabs((a).y))
    #define MAGMA_tally3_Z_CNJG(a)       cuConj(a)
    
    #define MAGMA_tally3_C_MAKE(r,i)     make_cuFloatComplex(r, i)
    #define MAGMA_tally3_C_REAL(a)       (a).x
    #define MAGMA_tally3_C_IMAG(a)       (a).y
    #define MAGMA_tally3_C_ADD(a, b)     cuCaddf(a, b)
    #define MAGMA_tally3_C_SUB(a, b)     cuCsubf(a, b)
    #define MAGMA_tally3_C_MUL(a, b)     cuCmulf(a, b)
    #define MAGMA_tally3_C_DIV(a, b)     cuCdivf(a, b)
    #define MAGMA_tally3_C_ABS(a)        cuCabsf(a)
    #define MAGMA_tally3_C_ABS1(a)       (fabsf((a).x) + fabsf((a).y))
    #define MAGMA_tally3_C_CNJG(a)       cuConjf(a)
    
#elif defined(HAVE_clAmdBlas)
    #include <clAmdBlas.h>
    
    typedef cl_command_queue  magma_tally3_queue_t;
    typedef cl_event          magma_tally3_event_t;
    typedef cl_device_id      magma_tally3_device_t;
    
    typedef DoubleComplex magma_tally3DoubleComplex;
    typedef FloatComplex  magma_tally3FloatComplex;
    
    #define MAGMA_tally3_Z_MAKE(r,i)     doubleComplex(r,i)
    #define MAGMA_tally3_Z_REAL(a)       (a).x
    #define MAGMA_tally3_Z_IMAG(a)       (a).y
    #define MAGMA_tally3_Z_ADD(a, b)     MAGMA_tally3_Z_MAKE((a).x+(b).x, (a).y+(b).y)
    #define MAGMA_tally3_Z_SUB(a, b)     MAGMA_tally3_Z_MAKE((a).x-(b).x, (a).y-(b).y)
    #define MAGMA_tally3_Z_ABS(a)        magma_tally3_cabs(a)
    #define MAGMA_tally3_Z_ABS1(a)       (fabs((a).x) + fabs((a).y))
    #define MAGMA_tally3_Z_CNJG(a)       MAGMA_tally3_Z_MAKE((a).x, -(a).y)
    
    #define MAGMA_tally3_C_MAKE(r,i)     floatComplex(r,i)
    #define MAGMA_tally3_C_REAL(a)       (a).x
    #define MAGMA_tally3_C_IMAG(a)       (a).y
    #define MAGMA_tally3_C_ADD(a, b)     MAGMA_tally3_C_MAKE((a).x+(b).x, (a).y+(b).y)
    #define MAGMA_tally3_C_SUB(a, b)     MAGMA_tally3_C_MAKE((a).x-(b).x, (a).y-(b).y)
    #define MAGMA_tally3_C_ABS(a)        magma_tally3_cabsf(a)
    #define MAGMA_tally3_C_ABS1(a)       (fabsf((a).x) + fabsf((a).y))
    #define MAGMA_tally3_C_CNJG(a)       MAGMA_tally3_C_MAKE((a).x, -(a).y)

#elif defined(HAVE_MIC)
    #include <stdio.h>
    #include <stdlib.h>
    #include <stdint.h>
    #include <unistd.h>
    #include <fcntl.h>
    #include <string.h>
    #include <sys/mman.h>
    #include <sys/ioctl.h>
    #include <sys/time.h>
    #include <scif.h>
    //#include <mkl.h>

    typedef int   magma_tally3_queue_t;
    typedef int   magma_tally3_event_t;
    typedef int   magma_tally3_device_t;

    #include <complex>
    typedef std::complex<float>   magma_tally3FloatComplex;
    typedef std::complex<double>  magma_tally3DoubleComplex;

    #define MAGMA_tally3_Z_MAKE(r, i)    std::complex<double>(r,i)
    #define MAGMA_tally3_Z_REAL(x)       (x).real()
    #define MAGMA_tally3_Z_IMAG(x)       (x).imag()
    #define MAGMA_tally3_Z_SET2REAL(a,r) { (a).real() = (r);   (a).imag() = 0.0; }
    #define MAGMA_tally3_Z_ADD(a, b)     ((a)+(b))
    #define MAGMA_tally3_Z_SUB(a, b)     ((a)-(b))
    #define MAGMA_tally3_Z_ABS(a)        abs(a)
    #define MAGMA_tally3_Z_MUL(a, b)     ((a)*(b))
    #define MAGMA_tally3_Z_DIV(a, b)     ((a)/(b))
    #define MAGMA_tally3_Z_ABS1(a)       (fabs((a).real()) + fabs((a).imag())) 
    #define MAGMA_tally3_Z_CNJG(a)       conj(a)
    #define MAGMA_tally3_Z_DSCALE(v,t,s) ((v) = (t)/(s))

    #define MAGMA_tally3_C_MAKE(r, i)    std::complex<float> (r,i)
    #define MAGMA_tally3_C_REAL(x)       (x).real()
    #define MAGMA_tally3_C_IMAG(x)       (x).imag()
    #define MAGMA_tally3_C_SET2REAL(a,r) { (a).real() = (r);   (a).imag() = 0.0; }
    #define MAGMA_tally3_C_ADD(a, b)     ((a)+(b))
    #define MAGMA_tally3_C_SUB(a, b)     ((a)-(b))
    #define MAGMA_tally3_C_ABS(a)        abs(a)
    #define MAGMA_tally3_C_MUL(a, b)     ((a)*(b))
    #define MAGMA_tally3_C_DIV(a, b)     ((a)/(b))
    #define MAGMA_tally3_C_ABS1(a)       (fabs((a).real()) + fabs((a).imag()))
    #define MAGMA_tally3_C_CNJG(a)       conj(a)
    #define MAGMA_tally3_C_SSCALE(v,t,s) ((v) = (t)/(s))
#else
    #error "One of HAVE_CUBLAS, HAVE_clAmdBlas, or HAVE_MIC must be defined. For example, add -DHAVE_CUBLAS to CFLAGS, or #define HAVE_CUBLAS before #include <magma_tally3.h>. In MAGMA_tally3, this happens in Makefile.internal."
#endif

#define MAGMA_tally3_Z_EQUAL(a,b)        (MAGMA_tally3_Z_REAL(a)==MAGMA_tally3_Z_REAL(b) && MAGMA_tally3_Z_IMAG(a)==MAGMA_tally3_Z_IMAG(b))
#define MAGMA_tally3_Z_NEGATE(a)         MAGMA_tally3_Z_MAKE( -MAGMA_tally3_Z_REAL(a), -MAGMA_tally3_Z_IMAG(a))

#define MAGMA_tally3_C_EQUAL(a,b)        (MAGMA_tally3_C_REAL(a)==MAGMA_tally3_C_REAL(b) && MAGMA_tally3_C_IMAG(a)==MAGMA_tally3_C_IMAG(b))
#define MAGMA_tally3_C_NEGATE(a)         MAGMA_tally3_C_MAKE( -MAGMA_tally3_C_REAL(a), -MAGMA_tally3_C_IMAG(a))

#define MAGMA_tally3_D_MAKE(r,i)         (r)
#define MAGMA_tally3_D_REAL(x)           (x)
#define MAGMA_tally3_D_IMAG(x)           (0.0)
#define MAGMA_tally3_D_SET2REAL(a,r)     (a) = (r)
#define MAGMA_tally3_D_ADD(a, b)         ((a) + (b))
#define MAGMA_tally3_D_SUB(a, b)         ((a) - (b))
#define MAGMA_tally3_D_MUL(a, b)         ((a) * (b))
#define MAGMA_tally3_D_DIV(a, b)         ((a) / (b))
#define MAGMA_tally3_D_ABS(a)            ((a)>0 ? (a) : -(a))
#define MAGMA_tally3_D_ABS1(a)           ((a)>0 ? (a) : -(a))
#define MAGMA_tally3_D_CNJG(a)           (a)
#define MAGMA_tally3_D_EQUAL(a,b)        ((a) == (b))
#define MAGMA_tally3_D_NEGATE(a)         (-a)
#define MAGMA_tally3_D_DSCALE(v, t, s)   (v) = (t)/(s)

#define MAGMA_tally3_S_MAKE(r,i)         (r)
#define MAGMA_tally3_S_REAL(x)           (x)
#define MAGMA_tally3_S_IMAG(x)           (0.0)
#define MAGMA_tally3_S_SET2REAL(a,r)     (a) = (r)
#define MAGMA_tally3_S_ADD(a, b)         ((a) + (b))
#define MAGMA_tally3_S_SUB(a, b)         ((a) - (b))
#define MAGMA_tally3_S_MUL(a, b)         ((a) * (b))
#define MAGMA_tally3_S_DIV(a, b)         ((a) / (b))
#define MAGMA_tally3_S_ABS(a)            ((a)>0 ? (a) : -(a))
#define MAGMA_tally3_S_ABS1(a)           ((a)>0 ? (a) : -(a))
#define MAGMA_tally3_S_CNJG(a)           (a)
#define MAGMA_tally3_S_EQUAL(a,b)        ((a) == (b))
#define MAGMA_tally3_S_NEGATE(a)         (-a)
#define MAGMA_tally3_S_SSCALE(v, t, s)   (v) = (t)/(s)

#define MAGMA_tally3_Z_ZERO              MAGMA_tally3_Z_MAKE( 0.0, 0.0)
#define MAGMA_tally3_Z_ONE               MAGMA_tally3_Z_MAKE( 1.0, 0.0)
#define MAGMA_tally3_Z_HALF              MAGMA_tally3_Z_MAKE( 0.5, 0.0)
#define MAGMA_tally3_Z_NEG_ONE           MAGMA_tally3_Z_MAKE(-1.0, 0.0)
#define MAGMA_tally3_Z_NEG_HALF          MAGMA_tally3_Z_MAKE(-0.5, 0.0)

#define MAGMA_tally3_C_ZERO              MAGMA_tally3_C_MAKE( 0.0, 0.0)
#define MAGMA_tally3_C_ONE               MAGMA_tally3_C_MAKE( 1.0, 0.0)
#define MAGMA_tally3_C_HALF              MAGMA_tally3_C_MAKE( 0.5, 0.0)
#define MAGMA_tally3_C_NEG_ONE           MAGMA_tally3_C_MAKE(-1.0, 0.0)
#define MAGMA_tally3_C_NEG_HALF          MAGMA_tally3_C_MAKE(-0.5, 0.0)

#define MAGMA_tally3_D_ZERO              ( 0.0)
#define MAGMA_tally3_D_ONE               ( 1.0)
#define MAGMA_tally3_D_HALF              ( 0.5)
#define MAGMA_tally3_D_NEG_ONE           (-1.0)
#define MAGMA_tally3_D_NEG_HALF          (-0.5)

#define MAGMA_tally3_S_ZERO              ( 0.0)
#define MAGMA_tally3_S_ONE               ( 1.0)
#define MAGMA_tally3_S_HALF              ( 0.5)
#define MAGMA_tally3_S_NEG_ONE           (-1.0)
#define MAGMA_tally3_S_NEG_HALF          (-0.5)

#ifndef CBLAS_SADDR
#define CBLAS_SADDR(a)  &(a)
#endif

#if defined(HAVE_clAmdBlas)
    // OpenCL uses opaque memory references on GPU
    typedef cl_mem magma_tally3_ptr;
    typedef cl_mem magma_tally3Int_ptr;
    typedef cl_mem magma_tally3Index_ptr;
    typedef cl_mem magma_tally3Float_ptr;
    typedef cl_mem magma_tally3Double_ptr;
    typedef cl_mem magma_tally3FloatComplex_ptr;
    typedef cl_mem magma_tally3DoubleComplex_ptr;
    
    typedef cl_mem magma_tally3_const_ptr;
    typedef cl_mem magma_tally3Int_const_ptr;
    typedef cl_mem magma_tally3Index_const_ptr;
    typedef cl_mem magma_tally3Float_const_ptr;
    typedef cl_mem magma_tally3Double_const_ptr;
    typedef cl_mem magma_tally3FloatComplex_const_ptr;
    typedef cl_mem magma_tally3DoubleComplex_const_ptr;
#else
    // MIC and CUDA use regular pointers on GPU
    typedef void               *magma_tally3_ptr;
    typedef magma_tally3_int_t        *magma_tally3Int_ptr;
    typedef magma_tally3_index_t      *magma_tally3Index_ptr;
    typedef float              *magma_tally3Float_ptr;
    typedef double             *magma_tally3Double_ptr;
    typedef magma_tally3FloatComplex  *magma_tally3FloatComplex_ptr;
    typedef magma_tally3DoubleComplex *magma_tally3DoubleComplex_ptr;
    
    typedef void               const *magma_tally3_const_ptr;
    typedef magma_tally3_int_t        const *magma_tally3Int_const_ptr;
    typedef magma_tally3_index_t      const *magma_tally3Index_const_ptr;
    typedef float              const *magma_tally3Float_const_ptr;
    typedef double             const *magma_tally3Double_const_ptr;
    typedef magma_tally3FloatComplex  const *magma_tally3FloatComplex_const_ptr;
    typedef magma_tally3DoubleComplex const *magma_tally3DoubleComplex_const_ptr;
#endif


// ========================================
// MAGMA_tally3 constants

// ----------------------------------------
#define MAGMA_tally3_VERSION_MAJOR 1
#define MAGMA_tally3_VERSION_MINOR 6
#define MAGMA_tally3_VERSION_MICRO 2

// stage is "svn", "beta#", "rc#" (release candidate), or blank ("") for final release
#define MAGMA_tally3_VERSION_STAGE ""

#define Magma_tally3MaxGPUs 8
#define Magma_tally3MaxDevices 8

// ----------------------------------------
// Return codes
// LAPACK argument errors are < 0 but > MAGMA_tally3_ERR.
// MAGMA_tally3 errors are < MAGMA_tally3_ERR.
#define MAGMA_tally3_SUCCESS               0
#define MAGMA_tally3_ERR                  -100
#define MAGMA_tally3_ERR_NOT_INITIALIZED  -101
#define MAGMA_tally3_ERR_REINITIALIZED    -102
#define MAGMA_tally3_ERR_NOT_SUPPORTED    -103
#define MAGMA_tally3_ERR_ILLEGAL_VALUE    -104
#define MAGMA_tally3_ERR_NOT_FOUND        -105
#define MAGMA_tally3_ERR_ALLOCATION       -106
#define MAGMA_tally3_ERR_INTERNAL_LIMIT   -107
#define MAGMA_tally3_ERR_UNALLOCATED      -108
#define MAGMA_tally3_ERR_FILESYSTEM       -109
#define MAGMA_tally3_ERR_UNEXPECTED       -110
#define MAGMA_tally3_ERR_SEQUENCE_FLUSHED -111
#define MAGMA_tally3_ERR_HOST_ALLOC       -112
#define MAGMA_tally3_ERR_DEVICE_ALLOC     -113
#define MAGMA_tally3_ERR_CUDASTREAM       -114
#define MAGMA_tally3_ERR_INVALID_PTR      -115
#define MAGMA_tally3_ERR_UNKNOWN          -116
#define MAGMA_tally3_ERR_NOT_IMPLEMENTED  -117

// some sparse-iter errors
#define MAGMA_tally3_SLOW_CONVERGENCE     -201
#define MAGMA_tally3_DIVERGENCE           -202
#define MAGMA_tally3_NONSPD               -203
#define MAGMA_tally3_ERR_BADPRECOND       -204

// When adding error codes, please add to interface_cuda/error.cpp

// map cusparse errors to magma_tally3 errors
#define MAGMA_tally3_ERR_CUSPARSE_NOT_INITIALIZED            -3001
#define MAGMA_tally3_ERR_CUSPARSE_ALLOC_FAILED               -3002
#define MAGMA_tally3_ERR_CUSPARSE_INVALID_VALUE              -3003
#define MAGMA_tally3_ERR_CUSPARSE_ARCH_MISMATCH              -3004
#define MAGMA_tally3_ERR_CUSPARSE_MAPPING_ERROR              -3005
#define MAGMA_tally3_ERR_CUSPARSE_EXECUTION_FAILED           -3006
#define MAGMA_tally3_ERR_CUSPARSE_INTERNAL_ERROR             -3007
#define MAGMA_tally3_ERR_CUSPARSE_MATRIX_TYPE_NOT_SUPPORTED  -3008
#define MAGMA_tally3_ERR_CUSPARSE_ZERO_PIVOT                 -3009


// ----------------------------------------
// parameter constants
// numbering is consistent with CBLAS and PLASMA; see plasma/include/plasma.h
// also with lapack_cwrapper/include/lapack_enum.h
typedef enum {
    Magma_tally3False         = 0,
    Magma_tally3True          = 1
} magma_tally3_bool_t;

typedef enum {
    Magma_tally3RowMajor      = 101,
    Magma_tally3ColMajor      = 102
} magma_tally3_order_t;

// Magma_tally3_ConjTrans is an alias for those rare occasions (zlarfb, zun*, zher*k)
// where we want Magma_tally3_ConjTrans to convert to Magma_tally3Trans in precision generation.
typedef enum {
    Magma_tally3NoTrans       = 111,
    Magma_tally3Trans         = 112,
    Magma_tally3ConjTrans     = 113,
    Magma_tally3_ConjTrans    = Magma_tally3ConjTrans
} magma_tally3_trans_t;

typedef enum {
    Magma_tally3Upper         = 121,
    Magma_tally3Lower         = 122,
    Magma_tally3UpperLower    = 123,
    Magma_tally3Full          = 123,  /* lascl, laset */
    Magma_tally3Hessenberg    = 124   /* lascl */
} magma_tally3_uplo_t;

typedef magma_tally3_uplo_t magma_tally3_type_t;  /* lascl */

typedef enum {
    Magma_tally3NonUnit       = 131,
    Magma_tally3Unit          = 132
} magma_tally3_diag_t;

typedef enum {
    Magma_tally3Left          = 141,
    Magma_tally3Right         = 142,
    Magma_tally3BothSides     = 143   /* trevc */
} magma_tally3_side_t;

typedef enum {
    Magma_tally3OneNorm       = 171,  /* lange, lanhe */
    Magma_tally3RealOneNorm   = 172,
    Magma_tally3TwoNorm       = 173,
    Magma_tally3FrobeniusNorm = 174,
    Magma_tally3InfNorm       = 175,
    Magma_tally3RealInfNorm   = 176,
    Magma_tally3MaxNorm       = 177,
    Magma_tally3RealMaxNorm   = 178
} magma_tally3_norm_t;

typedef enum {
    Magma_tally3DistUniform   = 201,  /* latms */
    Magma_tally3DistSymmetric = 202,
    Magma_tally3DistNormal    = 203
} magma_tally3_dist_t;

typedef enum {
    Magma_tally3HermGeev      = 241,  /* latms */
    Magma_tally3HermPoev      = 242,
    Magma_tally3NonsymPosv    = 243,
    Magma_tally3SymPosv       = 244
} magma_tally3_sym_t;

typedef enum {
    Magma_tally3NoPacking     = 291,  /* latms */
    Magma_tally3PackSubdiag   = 292,
    Magma_tally3PackSupdiag   = 293,
    Magma_tally3PackColumn    = 294,
    Magma_tally3PackRow       = 295,
    Magma_tally3PackLowerBand = 296,
    Magma_tally3PackUpeprBand = 297,
    Magma_tally3PackAll       = 298
} magma_tally3_pack_t;

typedef enum {
    Magma_tally3NoVec         = 301,  /* geev, syev, gesvd */
    Magma_tally3Vec           = 302,  /* geev, syev */
    Magma_tally3IVec          = 303,  /* stedc */
    Magma_tally3AllVec        = 304,  /* gesvd, trevc */
    Magma_tally3SomeVec       = 305,  /* gesvd, trevc */
    Magma_tally3OverwriteVec  = 306,  /* gesvd */
    Magma_tally3BacktransVec  = 307   /* trevc */
} magma_tally3_vec_t;

typedef enum {
    Magma_tally3RangeAll      = 311,  /* syevx, etc. */
    Magma_tally3RangeV        = 312,
    Magma_tally3RangeI        = 313
} magma_tally3_range_t;

typedef enum {
    Magma_tally3Q             = 322,  /* unmbr, ungbr */
    Magma_tally3P             = 323
} magma_tally3_vect_t;

typedef enum {
    Magma_tally3Forward       = 391,  /* larfb */
    Magma_tally3Backward      = 392
} magma_tally3_direct_t;

typedef enum {
    Magma_tally3Columnwise    = 401,  /* larfb */
    Magma_tally3Rowwise       = 402
} magma_tally3_storev_t;

// --------------------
// sparse
typedef enum {
    Magma_tally3_CSR          = 411,
    Magma_tally3_ELLPACKT     = 412,
    Magma_tally3_ELL          = 413,
    Magma_tally3_DENSE        = 414,
    Magma_tally3_BCSR         = 415,
    Magma_tally3_CSC          = 416,
    Magma_tally3_HYB          = 417,
    Magma_tally3_COO          = 418,
    Magma_tally3_ELLRT        = 419,
    Magma_tally3_SPMVFUNCTION = 420,
    Magma_tally3_SELLP        = 421,
    Magma_tally3_ELLD         = 422,

    Magma_tally3_CSRD         = 424,
    Magma_tally3_CSRL         = 427,
    Magma_tally3_CSRU         = 428,
    Magma_tally3_CSRCOO       = 429
} magma_tally3_storage_t;


typedef enum {
    Magma_tally3_CG           = 431,
    Magma_tally3_CGMERGE      = 432,
    Magma_tally3_GMRES        = 433,
    Magma_tally3_BICGSTAB     = 434,
  Magma_tally3_BICGSTABMERGE  = 435,
  Magma_tally3_BICGSTABMERGE2 = 436,
    Magma_tally3_JACOBI       = 437,
    Magma_tally3_GS           = 438,
    Magma_tally3_ITERREF      = 439,
    Magma_tally3_BCSRLU       = 440,
    Magma_tally3_PCG          = 441,
    Magma_tally3_PGMRES       = 442,
    Magma_tally3_PBICGSTAB    = 443,
    Magma_tally3_PASTIX       = 444,
    Magma_tally3_ILU          = 445,
    Magma_tally3_ICC          = 446,
    Magma_tally3_AILU         = 447,
    Magma_tally3_AICC         = 448,
    Magma_tally3_BAITER       = 449,
    Magma_tally3_LOBPCG       = 450,
    Magma_tally3_NONE         = 451,
    Magma_tally3_FUNCTION     = 452
} magma_tally3_solver_type;

typedef enum {
    Magma_tally3_CGS          = 461,
    Magma_tally3_FUSED_CGS    = 462,
    Magma_tally3_MGS          = 463
} magma_tally3_ortho_t;

typedef enum {
    Magma_tally3_CPU          = 471,
    Magma_tally3_DEV          = 472
} magma_tally3_location_t;

typedef enum {
    Magma_tally3_GENERAL      = 481,
    Magma_tally3_SYMMETRIC    = 482
} magma_tally3_symmetry_t;

typedef enum {
    Magma_tally3_ORDERED      = 491,
    Magma_tally3_DIAGFIRST    = 492,
    Magma_tally3_UNITY        = 493,
    Magma_tally3_VALUE        = 494
} magma_tally3_diagorder_t;

typedef enum {
    Magma_tally3_DCOMPLEX     = 501,
    Magma_tally3_FCOMPLEX     = 502,
    Magma_tally3_DOUBLE       = 503,
    Magma_tally3_FLOAT        = 504
} magma_tally3_precision;

typedef enum {
    Magma_tally3_NOSCALE      = 511,
    Magma_tally3_UNITROW      = 512,
    Magma_tally3_UNITDIAG     = 513
} magma_tally3_scale_t;

typedef enum {
    Magma_tally3_FULL         = 521,
    Magma_tally3_LOWER        = 522,
    Magma_tally3_UPPER        = 523
} magma_tally3_fillmode_t;



// When adding constants, remember to do these steps as appropriate:
// 1)  add magma_tally3_xxxx_const()  converter below and in control/constants.cpp
// 2a) add to magma_tally32lapack_const_tally3ants[] in control/constants.cpp
// 2b) update min & max here, which are used to check bounds for magma_tally32lapack_const_tally3ants[]
// 2c) add lapack_xxxx_const() converter below and in control/constants.cpp
#define Magma_tally32lapack_Min  Magma_tally3False     // 0
#define Magma_tally32lapack_Max  Magma_tally3Rowwise   // 402


// ----------------------------------------
// string constants for calling Fortran BLAS and LAPACK
// todo: use translators instead? lapack_const_tally3( Magma_tally3Upper )
#define Magma_tally3RowMajorStr      "Row"
#define Magma_tally3ColMajorStr      "Col"

#define Magma_tally3NoTransStr       "NoTrans"
#define Magma_tally3TransStr         "Trans"
#define Magma_tally3ConjTransStr     "ConjTrans"
#define Magma_tally3_ConjTransStr    "ConjTrans"

#define Magma_tally3UpperStr         "Upper"
#define Magma_tally3LowerStr         "Lower"
#define Magma_tally3UpperLowerStr    "Full"
#define Magma_tally3FullStr          "Full"

#define Magma_tally3NonUnitStr       "NonUnit"
#define Magma_tally3UnitStr          "Unit"

#define Magma_tally3LeftStr          "Left"
#define Magma_tally3RightStr         "Right"
#define Magma_tally3BothSidesStr     "Both"

#define Magma_tally3OneNormStr       "1"
#define Magma_tally3TwoNormStr       "2"
#define Magma_tally3FrobeniusNormStr "Fro"
#define Magma_tally3InfNormStr       "Inf"
#define Magma_tally3MaxNormStr       "Max"

#define Magma_tally3ForwardStr       "Forward"
#define Magma_tally3BackwardStr      "Backward"

#define Magma_tally3ColumnwiseStr    "Columnwise"
#define Magma_tally3RowwiseStr       "Rowwise"

#define Magma_tally3NoVecStr         "NoVec"
#define Magma_tally3VecStr           "Vec"
#define Magma_tally3IVecStr          "IVec"
#define Magma_tally3AllVecStr        "All"
#define Magma_tally3SomeVecStr       "Some"
#define Magma_tally3OverwriteVecStr  "Overwrite"


#ifdef __cplusplus
extern "C" {
#endif

// --------------------
// Convert LAPACK character constants to MAGMA_tally3 constants.
// This is a one-to-many mapping, requiring multiple translators
// (e.g., "N" can be NoTrans or NonUnit or NoVec).
magma_tally3_bool_t   magma_tally3_bool_const  ( char lapack_char );
magma_tally3_order_t  magma_tally3_order_const ( char lapack_char );
magma_tally3_trans_t  magma_tally3_trans_const ( char lapack_char );
magma_tally3_uplo_t   magma_tally3_uplo_const  ( char lapack_char );
magma_tally3_diag_t   magma_tally3_diag_const  ( char lapack_char );
magma_tally3_side_t   magma_tally3_side_const  ( char lapack_char );
magma_tally3_norm_t   magma_tally3_norm_const  ( char lapack_char );
magma_tally3_dist_t   magma_tally3_dist_const  ( char lapack_char );
magma_tally3_sym_t    magma_tally3_sym_const   ( char lapack_char );
magma_tally3_pack_t   magma_tally3_pack_const  ( char lapack_char );
magma_tally3_vec_t    magma_tally3_vec_const   ( char lapack_char );
magma_tally3_range_t  magma_tally3_range_const ( char lapack_char );
magma_tally3_vect_t   magma_tally3_vect_const  ( char lapack_char );
magma_tally3_direct_t magma_tally3_direct_const( char lapack_char );
magma_tally3_storev_t magma_tally3_storev_const( char lapack_char );


// --------------------
// Convert MAGMA_tally3 constants to LAPACK(E) constants.
// The generic lapack_const_tally3 works for all cases, but the specific routines
// (e.g., lapack_trans_const_tally3) do better error checking.
const char* lapack_const_tally3       ( int            magma_tally3_const );
const char* lapack_bool_const_tally3  ( magma_tally3_bool_t   magma_tally3_const );
const char* lapack_order_const_tally3 ( magma_tally3_order_t  magma_tally3_const );
const char* lapack_trans_const_tally3 ( magma_tally3_trans_t  magma_tally3_const );
const char* lapack_uplo_const_tally3  ( magma_tally3_uplo_t   magma_tally3_const );
const char* lapack_diag_const_tally3  ( magma_tally3_diag_t   magma_tally3_const );
const char* lapack_side_const_tally3  ( magma_tally3_side_t   magma_tally3_const );
const char* lapack_norm_const_tally3  ( magma_tally3_norm_t   magma_tally3_const );
const char* lapack_dist_const_tally3  ( magma_tally3_dist_t   magma_tally3_const );
const char* lapack_sym_const_tally3   ( magma_tally3_sym_t    magma_tally3_const );
const char* lapack_pack_const_tally3  ( magma_tally3_pack_t   magma_tally3_const );
const char* lapack_vec_const_tally3   ( magma_tally3_vec_t    magma_tally3_const );
const char* lapack_range_const_tally3 ( magma_tally3_range_t  magma_tally3_const );
const char* lapack_vect_const_tally3  ( magma_tally3_vect_t   magma_tally3_const );
const char* lapack_direct_const_tally3( magma_tally3_direct_t magma_tally3_const );
const char* lapack_storev_const_tally3( magma_tally3_storev_t magma_tally3_const );

static inline char lapacke_const_tally3       ( int magma_tally3_const            ) { return *lapack_const_tally3       ( magma_tally3_const ); }
static inline char lapacke_bool_const_tally3  ( magma_tally3_bool_t   magma_tally3_const ) { return *lapack_bool_const_tally3  ( magma_tally3_const ); }
static inline char lapacke_order_const_tally3 ( magma_tally3_order_t  magma_tally3_const ) { return *lapack_order_const_tally3 ( magma_tally3_const ); }
static inline char lapacke_trans_const_tally3 ( magma_tally3_trans_t  magma_tally3_const ) { return *lapack_trans_const_tally3 ( magma_tally3_const ); }
static inline char lapacke_uplo_const_tally3  ( magma_tally3_uplo_t   magma_tally3_const ) { return *lapack_uplo_const_tally3  ( magma_tally3_const ); }
static inline char lapacke_diag_const_tally3  ( magma_tally3_diag_t   magma_tally3_const ) { return *lapack_diag_const_tally3  ( magma_tally3_const ); }
static inline char lapacke_side_const_tally3  ( magma_tally3_side_t   magma_tally3_const ) { return *lapack_side_const_tally3  ( magma_tally3_const ); }
static inline char lapacke_norm_const_tally3  ( magma_tally3_norm_t   magma_tally3_const ) { return *lapack_norm_const_tally3  ( magma_tally3_const ); }
static inline char lapacke_dist_const_tally3  ( magma_tally3_dist_t   magma_tally3_const ) { return *lapack_dist_const_tally3  ( magma_tally3_const ); }
static inline char lapacke_sym_const_tally3   ( magma_tally3_sym_t    magma_tally3_const ) { return *lapack_sym_const_tally3   ( magma_tally3_const ); }
static inline char lapacke_pack_const_tally3  ( magma_tally3_pack_t   magma_tally3_const ) { return *lapack_pack_const_tally3  ( magma_tally3_const ); }
static inline char lapacke_vec_const_tally3   ( magma_tally3_vec_t    magma_tally3_const ) { return *lapack_vec_const_tally3   ( magma_tally3_const ); }
static inline char lapacke_range_const_tally3 ( magma_tally3_range_t  magma_tally3_const ) { return *lapack_range_const_tally3 ( magma_tally3_const ); }
static inline char lapacke_vect_const_tally3  ( magma_tally3_vect_t   magma_tally3_const ) { return *lapack_vect_const_tally3  ( magma_tally3_const ); }
static inline char lapacke_direct_const_tally3( magma_tally3_direct_t magma_tally3_const ) { return *lapack_direct_const_tally3( magma_tally3_const ); }
static inline char lapacke_storev_const_tally3( magma_tally3_storev_t magma_tally3_const ) { return *lapack_storev_const_tally3( magma_tally3_const ); }


// --------------------
// Convert MAGMA_tally3 constants to clAmdBlas constants.
#if defined(HAVE_clAmdBlas)
clAmdBlasOrder       amdblas_order_const( magma_tally3_order_t order );
clAmdBlasTranspose   amdblas_trans_const( magma_tally3_trans_t trans );
clAmdBlasUplo        amdblas_uplo_const ( magma_tally3_uplo_t  uplo  );
clAmdBlasDiag        amdblas_diag_const ( magma_tally3_diag_t  diag  );
clAmdBlasSide        amdblas_side_const ( magma_tally3_side_t  side  );
#endif


// --------------------
// Convert MAGMA_tally3 constants to CUBLAS constants.
#if defined(CUBLAS_V2_H_)
cublasOperation_t    cublas_trans_const_tally3 ( magma_tally3_trans_t trans );
cublasFillMode_t     cublas_uplo_const_tally3  ( magma_tally3_uplo_t  uplo  );
cublasDiagType_t     cublas_diag_const_tally3  ( magma_tally3_diag_t  diag  );
cublasSideMode_t     cublas_side_const_tally3  ( magma_tally3_side_t  side  );
#endif


// --------------------
// Convert MAGMA_tally3 constants to CBLAS constants.
#if defined(HAVE_CBLAS)
#include <cblas.h>
enum CBLAS_ORDER     cblas_order_const  ( magma_tally3_order_t order );
enum CBLAS_TRANSPOSE cblas_trans_const  ( magma_tally3_trans_t trans );
enum CBLAS_UPLO      cblas_uplo_const   ( magma_tally3_uplo_t  uplo  );
enum CBLAS_DIAG      cblas_diag_const   ( magma_tally3_diag_t  diag  );
enum CBLAS_SIDE      cblas_side_const   ( magma_tally3_side_t  side  );
#endif


#ifdef __cplusplus
}
#endif

#endif        //  #ifndef MAGMA_tally3_TYPES_H
