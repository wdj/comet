/*
    -- MAGMA_tally2 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015
*/

#ifndef MAGMA_tally2_TYPES_H
#define MAGMA_tally2_TYPES_H

#include <stdint.h>
#include <assert.h>


// each implementation of MAGMA_tally2 defines HAVE_* appropriately.
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
// Similar to magma_tally2_int_t we declare magma_tally2_index_t used for row/column indices in sparse
#if defined(MAGMA_tally2_ILP64) || defined(MKL_ILP64)
//typedef int64_t magma_tally2_int_t;
typedef long long int magma_tally2_int_t;  // MKL uses long long int, not int64_t
#else
typedef int magma_tally2_int_t;
#endif

typedef int magma_tally2_index_t;

// Define new type that the precision generator will not change (matches PLASMA)
typedef double real_Double_t;


// ========================================
// define types specific to implementation (CUDA, OpenCL, MIC)
// define macros to deal with complex numbers
#if defined(HAVE_CUBLAS)
    #ifndef CUBLAS_V2_H_
    #include <cublas.h>
    #endif
    
    typedef cudaStream_t   magma_tally2_queue_t;
    typedef cudaEvent_t    magma_tally2_event_t;
    typedef int            magma_tally2_device_t;
    
    typedef cuDoubleComplex magma_tally2DoubleComplex;
    typedef cuFloatComplex  magma_tally2FloatComplex;
    
    #define MAGMA_tally2_Z_MAKE(r,i)     make_cuDoubleComplex(r, i)
    #define MAGMA_tally2_Z_REAL(a)       (a).x
    #define MAGMA_tally2_Z_IMAG(a)       (a).y
    #define MAGMA_tally2_Z_ADD(a, b)     cuCadd(a, b)
    #define MAGMA_tally2_Z_SUB(a, b)     cuCsub(a, b)
    #define MAGMA_tally2_Z_MUL(a, b)     cuCmul(a, b)
    #define MAGMA_tally2_Z_DIV(a, b)     cuCdiv(a, b)
    #define MAGMA_tally2_Z_ABS(a)        cuCabs(a)
    #define MAGMA_tally2_Z_ABS1(a)       (fabs((a).x) + fabs((a).y))
    #define MAGMA_tally2_Z_CNJG(a)       cuConj(a)
    
    #define MAGMA_tally2_C_MAKE(r,i)     make_cuFloatComplex(r, i)
    #define MAGMA_tally2_C_REAL(a)       (a).x
    #define MAGMA_tally2_C_IMAG(a)       (a).y
    #define MAGMA_tally2_C_ADD(a, b)     cuCaddf(a, b)
    #define MAGMA_tally2_C_SUB(a, b)     cuCsubf(a, b)
    #define MAGMA_tally2_C_MUL(a, b)     cuCmulf(a, b)
    #define MAGMA_tally2_C_DIV(a, b)     cuCdivf(a, b)
    #define MAGMA_tally2_C_ABS(a)        cuCabsf(a)
    #define MAGMA_tally2_C_ABS1(a)       (fabsf((a).x) + fabsf((a).y))
    #define MAGMA_tally2_C_CNJG(a)       cuConjf(a)
    
#elif defined(HAVE_clAmdBlas)
    #include <clAmdBlas.h>
    
    typedef cl_command_queue  magma_tally2_queue_t;
    typedef cl_event          magma_tally2_event_t;
    typedef cl_device_id      magma_tally2_device_t;
    
    typedef DoubleComplex magma_tally2DoubleComplex;
    typedef FloatComplex  magma_tally2FloatComplex;
    
    #define MAGMA_tally2_Z_MAKE(r,i)     doubleComplex(r,i)
    #define MAGMA_tally2_Z_REAL(a)       (a).x
    #define MAGMA_tally2_Z_IMAG(a)       (a).y
    #define MAGMA_tally2_Z_ADD(a, b)     MAGMA_tally2_Z_MAKE((a).x+(b).x, (a).y+(b).y)
    #define MAGMA_tally2_Z_SUB(a, b)     MAGMA_tally2_Z_MAKE((a).x-(b).x, (a).y-(b).y)
    #define MAGMA_tally2_Z_ABS(a)        magma_tally2_cabs(a)
    #define MAGMA_tally2_Z_ABS1(a)       (fabs((a).x) + fabs((a).y))
    #define MAGMA_tally2_Z_CNJG(a)       MAGMA_tally2_Z_MAKE((a).x, -(a).y)
    
    #define MAGMA_tally2_C_MAKE(r,i)     floatComplex(r,i)
    #define MAGMA_tally2_C_REAL(a)       (a).x
    #define MAGMA_tally2_C_IMAG(a)       (a).y
    #define MAGMA_tally2_C_ADD(a, b)     MAGMA_tally2_C_MAKE((a).x+(b).x, (a).y+(b).y)
    #define MAGMA_tally2_C_SUB(a, b)     MAGMA_tally2_C_MAKE((a).x-(b).x, (a).y-(b).y)
    #define MAGMA_tally2_C_ABS(a)        magma_tally2_cabsf(a)
    #define MAGMA_tally2_C_ABS1(a)       (fabsf((a).x) + fabsf((a).y))
    #define MAGMA_tally2_C_CNJG(a)       MAGMA_tally2_C_MAKE((a).x, -(a).y)

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

    typedef int   magma_tally2_queue_t;
    typedef int   magma_tally2_event_t;
    typedef int   magma_tally2_device_t;

    #include <complex>
    typedef std::complex<float>   magma_tally2FloatComplex;
    typedef std::complex<double>  magma_tally2DoubleComplex;

    #define MAGMA_tally2_Z_MAKE(r, i)    std::complex<double>(r,i)
    #define MAGMA_tally2_Z_REAL(x)       (x).real()
    #define MAGMA_tally2_Z_IMAG(x)       (x).imag()
    #define MAGMA_tally2_Z_SET2REAL(a,r) { (a).real() = (r);   (a).imag() = 0.0; }
    #define MAGMA_tally2_Z_ADD(a, b)     ((a)+(b))
    #define MAGMA_tally2_Z_SUB(a, b)     ((a)-(b))
    #define MAGMA_tally2_Z_ABS(a)        abs(a)
    #define MAGMA_tally2_Z_MUL(a, b)     ((a)*(b))
    #define MAGMA_tally2_Z_DIV(a, b)     ((a)/(b))
    #define MAGMA_tally2_Z_ABS1(a)       (fabs((a).real()) + fabs((a).imag())) 
    #define MAGMA_tally2_Z_CNJG(a)       conj(a)
    #define MAGMA_tally2_Z_DSCALE(v,t,s) ((v) = (t)/(s))

    #define MAGMA_tally2_C_MAKE(r, i)    std::complex<float> (r,i)
    #define MAGMA_tally2_C_REAL(x)       (x).real()
    #define MAGMA_tally2_C_IMAG(x)       (x).imag()
    #define MAGMA_tally2_C_SET2REAL(a,r) { (a).real() = (r);   (a).imag() = 0.0; }
    #define MAGMA_tally2_C_ADD(a, b)     ((a)+(b))
    #define MAGMA_tally2_C_SUB(a, b)     ((a)-(b))
    #define MAGMA_tally2_C_ABS(a)        abs(a)
    #define MAGMA_tally2_C_MUL(a, b)     ((a)*(b))
    #define MAGMA_tally2_C_DIV(a, b)     ((a)/(b))
    #define MAGMA_tally2_C_ABS1(a)       (fabs((a).real()) + fabs((a).imag()))
    #define MAGMA_tally2_C_CNJG(a)       conj(a)
    #define MAGMA_tally2_C_SSCALE(v,t,s) ((v) = (t)/(s))
#else
    #error "One of HAVE_CUBLAS, HAVE_clAmdBlas, or HAVE_MIC must be defined. For example, add -DHAVE_CUBLAS to CFLAGS, or #define HAVE_CUBLAS before #include <magma_tally2.h>. In MAGMA_tally2, this happens in Makefile.internal."
#endif

#define MAGMA_tally2_Z_EQUAL(a,b)        (MAGMA_tally2_Z_REAL(a)==MAGMA_tally2_Z_REAL(b) && MAGMA_tally2_Z_IMAG(a)==MAGMA_tally2_Z_IMAG(b))
#define MAGMA_tally2_Z_NEGATE(a)         MAGMA_tally2_Z_MAKE( -MAGMA_tally2_Z_REAL(a), -MAGMA_tally2_Z_IMAG(a))

#define MAGMA_tally2_C_EQUAL(a,b)        (MAGMA_tally2_C_REAL(a)==MAGMA_tally2_C_REAL(b) && MAGMA_tally2_C_IMAG(a)==MAGMA_tally2_C_IMAG(b))
#define MAGMA_tally2_C_NEGATE(a)         MAGMA_tally2_C_MAKE( -MAGMA_tally2_C_REAL(a), -MAGMA_tally2_C_IMAG(a))

#define MAGMA_tally2_D_MAKE(r,i)         (r)
#define MAGMA_tally2_D_REAL(x)           (x)
#define MAGMA_tally2_D_IMAG(x)           (0.0)
#define MAGMA_tally2_D_SET2REAL(a,r)     (a) = (r)
#define MAGMA_tally2_D_ADD(a, b)         ((a) + (b))
#define MAGMA_tally2_D_SUB(a, b)         ((a) - (b))
#define MAGMA_tally2_D_MUL(a, b)         ((a) * (b))
#define MAGMA_tally2_D_DIV(a, b)         ((a) / (b))
#define MAGMA_tally2_D_ABS(a)            ((a)>0 ? (a) : -(a))
#define MAGMA_tally2_D_ABS1(a)           ((a)>0 ? (a) : -(a))
#define MAGMA_tally2_D_CNJG(a)           (a)
#define MAGMA_tally2_D_EQUAL(a,b)        ((a) == (b))
#define MAGMA_tally2_D_NEGATE(a)         (-a)
#define MAGMA_tally2_D_DSCALE(v, t, s)   (v) = (t)/(s)

#define MAGMA_tally2_S_MAKE(r,i)         (r)
#define MAGMA_tally2_S_REAL(x)           (x)
#define MAGMA_tally2_S_IMAG(x)           (0.0)
#define MAGMA_tally2_S_SET2REAL(a,r)     (a) = (r)
#define MAGMA_tally2_S_ADD(a, b)         ((a) + (b))
#define MAGMA_tally2_S_SUB(a, b)         ((a) - (b))
#define MAGMA_tally2_S_MUL(a, b)         ((a) * (b))
#define MAGMA_tally2_S_DIV(a, b)         ((a) / (b))
#define MAGMA_tally2_S_ABS(a)            ((a)>0 ? (a) : -(a))
#define MAGMA_tally2_S_ABS1(a)           ((a)>0 ? (a) : -(a))
#define MAGMA_tally2_S_CNJG(a)           (a)
#define MAGMA_tally2_S_EQUAL(a,b)        ((a) == (b))
#define MAGMA_tally2_S_NEGATE(a)         (-a)
#define MAGMA_tally2_S_SSCALE(v, t, s)   (v) = (t)/(s)

#define MAGMA_tally2_Z_ZERO              MAGMA_tally2_Z_MAKE( 0.0, 0.0)
#define MAGMA_tally2_Z_ONE               MAGMA_tally2_Z_MAKE( 1.0, 0.0)
#define MAGMA_tally2_Z_HALF              MAGMA_tally2_Z_MAKE( 0.5, 0.0)
#define MAGMA_tally2_Z_NEG_ONE           MAGMA_tally2_Z_MAKE(-1.0, 0.0)
#define MAGMA_tally2_Z_NEG_HALF          MAGMA_tally2_Z_MAKE(-0.5, 0.0)

#define MAGMA_tally2_C_ZERO              MAGMA_tally2_C_MAKE( 0.0, 0.0)
#define MAGMA_tally2_C_ONE               MAGMA_tally2_C_MAKE( 1.0, 0.0)
#define MAGMA_tally2_C_HALF              MAGMA_tally2_C_MAKE( 0.5, 0.0)
#define MAGMA_tally2_C_NEG_ONE           MAGMA_tally2_C_MAKE(-1.0, 0.0)
#define MAGMA_tally2_C_NEG_HALF          MAGMA_tally2_C_MAKE(-0.5, 0.0)

#define MAGMA_tally2_D_ZERO              ( 0.0)
#define MAGMA_tally2_D_ONE               ( 1.0)
#define MAGMA_tally2_D_HALF              ( 0.5)
#define MAGMA_tally2_D_NEG_ONE           (-1.0)
#define MAGMA_tally2_D_NEG_HALF          (-0.5)

#define MAGMA_tally2_S_ZERO              ( 0.0)
#define MAGMA_tally2_S_ONE               ( 1.0)
#define MAGMA_tally2_S_HALF              ( 0.5)
#define MAGMA_tally2_S_NEG_ONE           (-1.0)
#define MAGMA_tally2_S_NEG_HALF          (-0.5)

#ifndef CBLAS_SADDR
#define CBLAS_SADDR(a)  &(a)
#endif

#if defined(HAVE_clAmdBlas)
    // OpenCL uses opaque memory references on GPU
    typedef cl_mem magma_tally2_ptr;
    typedef cl_mem magma_tally2Int_ptr;
    typedef cl_mem magma_tally2Index_ptr;
    typedef cl_mem magma_tally2Float_ptr;
    typedef cl_mem magma_tally2Double_ptr;
    typedef cl_mem magma_tally2FloatComplex_ptr;
    typedef cl_mem magma_tally2DoubleComplex_ptr;
    
    typedef cl_mem magma_tally2_const_ptr;
    typedef cl_mem magma_tally2Int_const_ptr;
    typedef cl_mem magma_tally2Index_const_ptr;
    typedef cl_mem magma_tally2Float_const_ptr;
    typedef cl_mem magma_tally2Double_const_ptr;
    typedef cl_mem magma_tally2FloatComplex_const_ptr;
    typedef cl_mem magma_tally2DoubleComplex_const_ptr;
#else
    // MIC and CUDA use regular pointers on GPU
    typedef void               *magma_tally2_ptr;
    typedef magma_tally2_int_t        *magma_tally2Int_ptr;
    typedef magma_tally2_index_t      *magma_tally2Index_ptr;
    typedef float              *magma_tally2Float_ptr;
    typedef double             *magma_tally2Double_ptr;
    typedef magma_tally2FloatComplex  *magma_tally2FloatComplex_ptr;
    typedef magma_tally2DoubleComplex *magma_tally2DoubleComplex_ptr;
    
    typedef void               const *magma_tally2_const_ptr;
    typedef magma_tally2_int_t        const *magma_tally2Int_const_ptr;
    typedef magma_tally2_index_t      const *magma_tally2Index_const_ptr;
    typedef float              const *magma_tally2Float_const_ptr;
    typedef double             const *magma_tally2Double_const_ptr;
    typedef magma_tally2FloatComplex  const *magma_tally2FloatComplex_const_ptr;
    typedef magma_tally2DoubleComplex const *magma_tally2DoubleComplex_const_ptr;
#endif


// ========================================
// MAGMA_tally2 constants

// ----------------------------------------
#define MAGMA_tally2_VERSION_MAJOR 1
#define MAGMA_tally2_VERSION_MINOR 6
#define MAGMA_tally2_VERSION_MICRO 2

// stage is "svn", "beta#", "rc#" (release candidate), or blank ("") for final release
#define MAGMA_tally2_VERSION_STAGE ""

#define Magma_tally2MaxGPUs 8
#define Magma_tally2MaxDevices 8

// ----------------------------------------
// Return codes
// LAPACK argument errors are < 0 but > MAGMA_tally2_ERR.
// MAGMA_tally2 errors are < MAGMA_tally2_ERR.
#define MAGMA_tally2_SUCCESS               0
#define MAGMA_tally2_ERR                  -100
#define MAGMA_tally2_ERR_NOT_INITIALIZED  -101
#define MAGMA_tally2_ERR_REINITIALIZED    -102
#define MAGMA_tally2_ERR_NOT_SUPPORTED    -103
#define MAGMA_tally2_ERR_ILLEGAL_VALUE    -104
#define MAGMA_tally2_ERR_NOT_FOUND        -105
#define MAGMA_tally2_ERR_ALLOCATION       -106
#define MAGMA_tally2_ERR_INTERNAL_LIMIT   -107
#define MAGMA_tally2_ERR_UNALLOCATED      -108
#define MAGMA_tally2_ERR_FILESYSTEM       -109
#define MAGMA_tally2_ERR_UNEXPECTED       -110
#define MAGMA_tally2_ERR_SEQUENCE_FLUSHED -111
#define MAGMA_tally2_ERR_HOST_ALLOC       -112
#define MAGMA_tally2_ERR_DEVICE_ALLOC     -113
#define MAGMA_tally2_ERR_CUDASTREAM       -114
#define MAGMA_tally2_ERR_INVALID_PTR      -115
#define MAGMA_tally2_ERR_UNKNOWN          -116
#define MAGMA_tally2_ERR_NOT_IMPLEMENTED  -117

// some sparse-iter errors
#define MAGMA_tally2_SLOW_CONVERGENCE     -201
#define MAGMA_tally2_DIVERGENCE           -202
#define MAGMA_tally2_NONSPD               -203
#define MAGMA_tally2_ERR_BADPRECOND       -204

// When adding error codes, please add to interface_cuda/error.cpp

// map cusparse errors to magma_tally2 errors
#define MAGMA_tally2_ERR_CUSPARSE_NOT_INITIALIZED            -3001
#define MAGMA_tally2_ERR_CUSPARSE_ALLOC_FAILED               -3002
#define MAGMA_tally2_ERR_CUSPARSE_INVALID_VALUE              -3003
#define MAGMA_tally2_ERR_CUSPARSE_ARCH_MISMATCH              -3004
#define MAGMA_tally2_ERR_CUSPARSE_MAPPING_ERROR              -3005
#define MAGMA_tally2_ERR_CUSPARSE_EXECUTION_FAILED           -3006
#define MAGMA_tally2_ERR_CUSPARSE_INTERNAL_ERROR             -3007
#define MAGMA_tally2_ERR_CUSPARSE_MATRIX_TYPE_NOT_SUPPORTED  -3008
#define MAGMA_tally2_ERR_CUSPARSE_ZERO_PIVOT                 -3009


// ----------------------------------------
// parameter constants
// numbering is consistent with CBLAS and PLASMA; see plasma/include/plasma.h
// also with lapack_cwrapper/include/lapack_enum.h
typedef enum {
    Magma_tally2False         = 0,
    Magma_tally2True          = 1
} magma_tally2_bool_t;

typedef enum {
    Magma_tally2RowMajor      = 101,
    Magma_tally2ColMajor      = 102
} magma_tally2_order_t;

// Magma_tally2_ConjTrans is an alias for those rare occasions (zlarfb, zun*, zher*k)
// where we want Magma_tally2_ConjTrans to convert to Magma_tally2Trans in precision generation.
typedef enum {
    Magma_tally2NoTrans       = 111,
    Magma_tally2Trans         = 112,
    Magma_tally2ConjTrans     = 113,
    Magma_tally2_ConjTrans    = Magma_tally2ConjTrans
} magma_tally2_trans_t;

typedef enum {
    Magma_tally2Upper         = 121,
    Magma_tally2Lower         = 122,
    Magma_tally2UpperLower    = 123,
    Magma_tally2Full          = 123,  /* lascl, laset */
    Magma_tally2Hessenberg    = 124   /* lascl */
} magma_tally2_uplo_t;

typedef magma_tally2_uplo_t magma_tally2_type_t;  /* lascl */

typedef enum {
    Magma_tally2NonUnit       = 131,
    Magma_tally2Unit          = 132
} magma_tally2_diag_t;

typedef enum {
    Magma_tally2Left          = 141,
    Magma_tally2Right         = 142,
    Magma_tally2BothSides     = 143   /* trevc */
} magma_tally2_side_t;

typedef enum {
    Magma_tally2OneNorm       = 171,  /* lange, lanhe */
    Magma_tally2RealOneNorm   = 172,
    Magma_tally2TwoNorm       = 173,
    Magma_tally2FrobeniusNorm = 174,
    Magma_tally2InfNorm       = 175,
    Magma_tally2RealInfNorm   = 176,
    Magma_tally2MaxNorm       = 177,
    Magma_tally2RealMaxNorm   = 178
} magma_tally2_norm_t;

typedef enum {
    Magma_tally2DistUniform   = 201,  /* latms */
    Magma_tally2DistSymmetric = 202,
    Magma_tally2DistNormal    = 203
} magma_tally2_dist_t;

typedef enum {
    Magma_tally2HermGeev      = 241,  /* latms */
    Magma_tally2HermPoev      = 242,
    Magma_tally2NonsymPosv    = 243,
    Magma_tally2SymPosv       = 244
} magma_tally2_sym_t;

typedef enum {
    Magma_tally2NoPacking     = 291,  /* latms */
    Magma_tally2PackSubdiag   = 292,
    Magma_tally2PackSupdiag   = 293,
    Magma_tally2PackColumn    = 294,
    Magma_tally2PackRow       = 295,
    Magma_tally2PackLowerBand = 296,
    Magma_tally2PackUpeprBand = 297,
    Magma_tally2PackAll       = 298
} magma_tally2_pack_t;

typedef enum {
    Magma_tally2NoVec         = 301,  /* geev, syev, gesvd */
    Magma_tally2Vec           = 302,  /* geev, syev */
    Magma_tally2IVec          = 303,  /* stedc */
    Magma_tally2AllVec        = 304,  /* gesvd, trevc */
    Magma_tally2SomeVec       = 305,  /* gesvd, trevc */
    Magma_tally2OverwriteVec  = 306,  /* gesvd */
    Magma_tally2BacktransVec  = 307   /* trevc */
} magma_tally2_vec_t;

typedef enum {
    Magma_tally2RangeAll      = 311,  /* syevx, etc. */
    Magma_tally2RangeV        = 312,
    Magma_tally2RangeI        = 313
} magma_tally2_range_t;

typedef enum {
    Magma_tally2Q             = 322,  /* unmbr, ungbr */
    Magma_tally2P             = 323
} magma_tally2_vect_t;

typedef enum {
    Magma_tally2Forward       = 391,  /* larfb */
    Magma_tally2Backward      = 392
} magma_tally2_direct_t;

typedef enum {
    Magma_tally2Columnwise    = 401,  /* larfb */
    Magma_tally2Rowwise       = 402
} magma_tally2_storev_t;

// --------------------
// sparse
typedef enum {
    Magma_tally2_CSR          = 411,
    Magma_tally2_ELLPACKT     = 412,
    Magma_tally2_ELL          = 413,
    Magma_tally2_DENSE        = 414,
    Magma_tally2_BCSR         = 415,
    Magma_tally2_CSC          = 416,
    Magma_tally2_HYB          = 417,
    Magma_tally2_COO          = 418,
    Magma_tally2_ELLRT        = 419,
    Magma_tally2_SPMVFUNCTION = 420,
    Magma_tally2_SELLP        = 421,
    Magma_tally2_ELLD         = 422,

    Magma_tally2_CSRD         = 424,
    Magma_tally2_CSRL         = 427,
    Magma_tally2_CSRU         = 428,
    Magma_tally2_CSRCOO       = 429
} magma_tally2_storage_t;


typedef enum {
    Magma_tally2_CG           = 431,
    Magma_tally2_CGMERGE      = 432,
    Magma_tally2_GMRES        = 433,
    Magma_tally2_BICGSTAB     = 434,
  Magma_tally2_BICGSTABMERGE  = 435,
  Magma_tally2_BICGSTABMERGE2 = 436,
    Magma_tally2_JACOBI       = 437,
    Magma_tally2_GS           = 438,
    Magma_tally2_ITERREF      = 439,
    Magma_tally2_BCSRLU       = 440,
    Magma_tally2_PCG          = 441,
    Magma_tally2_PGMRES       = 442,
    Magma_tally2_PBICGSTAB    = 443,
    Magma_tally2_PASTIX       = 444,
    Magma_tally2_ILU          = 445,
    Magma_tally2_ICC          = 446,
    Magma_tally2_AILU         = 447,
    Magma_tally2_AICC         = 448,
    Magma_tally2_BAITER       = 449,
    Magma_tally2_LOBPCG       = 450,
    Magma_tally2_NONE         = 451,
    Magma_tally2_FUNCTION     = 452
} magma_tally2_solver_type;

typedef enum {
    Magma_tally2_CGS          = 461,
    Magma_tally2_FUSED_CGS    = 462,
    Magma_tally2_MGS          = 463
} magma_tally2_ortho_t;

typedef enum {
    Magma_tally2_CPU          = 471,
    Magma_tally2_DEV          = 472
} magma_tally2_location_t;

typedef enum {
    Magma_tally2_GENERAL      = 481,
    Magma_tally2_SYMMETRIC    = 482
} magma_tally2_symmetry_t;

typedef enum {
    Magma_tally2_ORDERED      = 491,
    Magma_tally2_DIAGFIRST    = 492,
    Magma_tally2_UNITY        = 493,
    Magma_tally2_VALUE        = 494
} magma_tally2_diagorder_t;

typedef enum {
    Magma_tally2_DCOMPLEX     = 501,
    Magma_tally2_FCOMPLEX     = 502,
    Magma_tally2_DOUBLE       = 503,
    Magma_tally2_FLOAT        = 504
} magma_tally2_precision;

typedef enum {
    Magma_tally2_NOSCALE      = 511,
    Magma_tally2_UNITROW      = 512,
    Magma_tally2_UNITDIAG     = 513
} magma_tally2_scale_t;

typedef enum {
    Magma_tally2_FULL         = 521,
    Magma_tally2_LOWER        = 522,
    Magma_tally2_UPPER        = 523
} magma_tally2_fillmode_t;



// When adding constants, remember to do these steps as appropriate:
// 1)  add magma_tally2_xxxx_const()  converter below and in control/constants.cpp
// 2a) add to magma_tally22lapack_const_tally2ants[] in control/constants.cpp
// 2b) update min & max here, which are used to check bounds for magma_tally22lapack_const_tally2ants[]
// 2c) add lapack_xxxx_const() converter below and in control/constants.cpp
#define Magma_tally22lapack_Min  Magma_tally2False     // 0
#define Magma_tally22lapack_Max  Magma_tally2Rowwise   // 402


// ----------------------------------------
// string constants for calling Fortran BLAS and LAPACK
// todo: use translators instead? lapack_const_tally2( Magma_tally2Upper )
#define Magma_tally2RowMajorStr      "Row"
#define Magma_tally2ColMajorStr      "Col"

#define Magma_tally2NoTransStr       "NoTrans"
#define Magma_tally2TransStr         "Trans"
#define Magma_tally2ConjTransStr     "ConjTrans"
#define Magma_tally2_ConjTransStr    "ConjTrans"

#define Magma_tally2UpperStr         "Upper"
#define Magma_tally2LowerStr         "Lower"
#define Magma_tally2UpperLowerStr    "Full"
#define Magma_tally2FullStr          "Full"

#define Magma_tally2NonUnitStr       "NonUnit"
#define Magma_tally2UnitStr          "Unit"

#define Magma_tally2LeftStr          "Left"
#define Magma_tally2RightStr         "Right"
#define Magma_tally2BothSidesStr     "Both"

#define Magma_tally2OneNormStr       "1"
#define Magma_tally2TwoNormStr       "2"
#define Magma_tally2FrobeniusNormStr "Fro"
#define Magma_tally2InfNormStr       "Inf"
#define Magma_tally2MaxNormStr       "Max"

#define Magma_tally2ForwardStr       "Forward"
#define Magma_tally2BackwardStr      "Backward"

#define Magma_tally2ColumnwiseStr    "Columnwise"
#define Magma_tally2RowwiseStr       "Rowwise"

#define Magma_tally2NoVecStr         "NoVec"
#define Magma_tally2VecStr           "Vec"
#define Magma_tally2IVecStr          "IVec"
#define Magma_tally2AllVecStr        "All"
#define Magma_tally2SomeVecStr       "Some"
#define Magma_tally2OverwriteVecStr  "Overwrite"


#ifdef __cplusplus
extern "C" {
#endif

// --------------------
// Convert LAPACK character constants to MAGMA_tally2 constants.
// This is a one-to-many mapping, requiring multiple translators
// (e.g., "N" can be NoTrans or NonUnit or NoVec).
magma_tally2_bool_t   magma_tally2_bool_const  ( char lapack_char );
magma_tally2_order_t  magma_tally2_order_const ( char lapack_char );
magma_tally2_trans_t  magma_tally2_trans_const ( char lapack_char );
magma_tally2_uplo_t   magma_tally2_uplo_const  ( char lapack_char );
magma_tally2_diag_t   magma_tally2_diag_const  ( char lapack_char );
magma_tally2_side_t   magma_tally2_side_const  ( char lapack_char );
magma_tally2_norm_t   magma_tally2_norm_const  ( char lapack_char );
magma_tally2_dist_t   magma_tally2_dist_const  ( char lapack_char );
magma_tally2_sym_t    magma_tally2_sym_const   ( char lapack_char );
magma_tally2_pack_t   magma_tally2_pack_const  ( char lapack_char );
magma_tally2_vec_t    magma_tally2_vec_const   ( char lapack_char );
magma_tally2_range_t  magma_tally2_range_const ( char lapack_char );
magma_tally2_vect_t   magma_tally2_vect_const  ( char lapack_char );
magma_tally2_direct_t magma_tally2_direct_const( char lapack_char );
magma_tally2_storev_t magma_tally2_storev_const( char lapack_char );


// --------------------
// Convert MAGMA_tally2 constants to LAPACK(E) constants.
// The generic lapack_const_tally2 works for all cases, but the specific routines
// (e.g., lapack_trans_const_tally2) do better error checking.
const char* lapack_const_tally2       ( int            magma_tally2_const );
const char* lapack_bool_const_tally2  ( magma_tally2_bool_t   magma_tally2_const );
const char* lapack_order_const_tally2 ( magma_tally2_order_t  magma_tally2_const );
const char* lapack_trans_const_tally2 ( magma_tally2_trans_t  magma_tally2_const );
const char* lapack_uplo_const_tally2  ( magma_tally2_uplo_t   magma_tally2_const );
const char* lapack_diag_const_tally2  ( magma_tally2_diag_t   magma_tally2_const );
const char* lapack_side_const_tally2  ( magma_tally2_side_t   magma_tally2_const );
const char* lapack_norm_const_tally2  ( magma_tally2_norm_t   magma_tally2_const );
const char* lapack_dist_const_tally2  ( magma_tally2_dist_t   magma_tally2_const );
const char* lapack_sym_const_tally2   ( magma_tally2_sym_t    magma_tally2_const );
const char* lapack_pack_const_tally2  ( magma_tally2_pack_t   magma_tally2_const );
const char* lapack_vec_const_tally2   ( magma_tally2_vec_t    magma_tally2_const );
const char* lapack_range_const_tally2 ( magma_tally2_range_t  magma_tally2_const );
const char* lapack_vect_const_tally2  ( magma_tally2_vect_t   magma_tally2_const );
const char* lapack_direct_const_tally2( magma_tally2_direct_t magma_tally2_const );
const char* lapack_storev_const_tally2( magma_tally2_storev_t magma_tally2_const );

static inline char lapacke_const_tally2       ( int magma_tally2_const            ) { return *lapack_const_tally2       ( magma_tally2_const ); }
static inline char lapacke_bool_const_tally2  ( magma_tally2_bool_t   magma_tally2_const ) { return *lapack_bool_const_tally2  ( magma_tally2_const ); }
static inline char lapacke_order_const_tally2 ( magma_tally2_order_t  magma_tally2_const ) { return *lapack_order_const_tally2 ( magma_tally2_const ); }
static inline char lapacke_trans_const_tally2 ( magma_tally2_trans_t  magma_tally2_const ) { return *lapack_trans_const_tally2 ( magma_tally2_const ); }
static inline char lapacke_uplo_const_tally2  ( magma_tally2_uplo_t   magma_tally2_const ) { return *lapack_uplo_const_tally2  ( magma_tally2_const ); }
static inline char lapacke_diag_const_tally2  ( magma_tally2_diag_t   magma_tally2_const ) { return *lapack_diag_const_tally2  ( magma_tally2_const ); }
static inline char lapacke_side_const_tally2  ( magma_tally2_side_t   magma_tally2_const ) { return *lapack_side_const_tally2  ( magma_tally2_const ); }
static inline char lapacke_norm_const_tally2  ( magma_tally2_norm_t   magma_tally2_const ) { return *lapack_norm_const_tally2  ( magma_tally2_const ); }
static inline char lapacke_dist_const_tally2  ( magma_tally2_dist_t   magma_tally2_const ) { return *lapack_dist_const_tally2  ( magma_tally2_const ); }
static inline char lapacke_sym_const_tally2   ( magma_tally2_sym_t    magma_tally2_const ) { return *lapack_sym_const_tally2   ( magma_tally2_const ); }
static inline char lapacke_pack_const_tally2  ( magma_tally2_pack_t   magma_tally2_const ) { return *lapack_pack_const_tally2  ( magma_tally2_const ); }
static inline char lapacke_vec_const_tally2   ( magma_tally2_vec_t    magma_tally2_const ) { return *lapack_vec_const_tally2   ( magma_tally2_const ); }
static inline char lapacke_range_const_tally2 ( magma_tally2_range_t  magma_tally2_const ) { return *lapack_range_const_tally2 ( magma_tally2_const ); }
static inline char lapacke_vect_const_tally2  ( magma_tally2_vect_t   magma_tally2_const ) { return *lapack_vect_const_tally2  ( magma_tally2_const ); }
static inline char lapacke_direct_const_tally2( magma_tally2_direct_t magma_tally2_const ) { return *lapack_direct_const_tally2( magma_tally2_const ); }
static inline char lapacke_storev_const_tally2( magma_tally2_storev_t magma_tally2_const ) { return *lapack_storev_const_tally2( magma_tally2_const ); }


// --------------------
// Convert MAGMA_tally2 constants to clAmdBlas constants.
#if defined(HAVE_clAmdBlas)
clAmdBlasOrder       amdblas_order_const( magma_tally2_order_t order );
clAmdBlasTranspose   amdblas_trans_const( magma_tally2_trans_t trans );
clAmdBlasUplo        amdblas_uplo_const ( magma_tally2_uplo_t  uplo  );
clAmdBlasDiag        amdblas_diag_const ( magma_tally2_diag_t  diag  );
clAmdBlasSide        amdblas_side_const ( magma_tally2_side_t  side  );
#endif


// --------------------
// Convert MAGMA_tally2 constants to CUBLAS constants.
#if defined(CUBLAS_V2_H_)
cublasOperation_t    cublas_trans_const_tally2 ( magma_tally2_trans_t trans );
cublasFillMode_t     cublas_uplo_const_tally2  ( magma_tally2_uplo_t  uplo  );
cublasDiagType_t     cublas_diag_const_tally2  ( magma_tally2_diag_t  diag  );
cublasSideMode_t     cublas_side_const_tally2  ( magma_tally2_side_t  side  );
#endif


// --------------------
// Convert MAGMA_tally2 constants to CBLAS constants.
#if defined(HAVE_CBLAS)
#include <cblas.h>
enum CBLAS_ORDER     cblas_order_const  ( magma_tally2_order_t order );
enum CBLAS_TRANSPOSE cblas_trans_const  ( magma_tally2_trans_t trans );
enum CBLAS_UPLO      cblas_uplo_const   ( magma_tally2_uplo_t  uplo  );
enum CBLAS_DIAG      cblas_diag_const   ( magma_tally2_diag_t  diag  );
enum CBLAS_SIDE      cblas_side_const   ( magma_tally2_side_t  side  );
#endif


#ifdef __cplusplus
}
#endif

#endif        //  #ifndef MAGMA_tally2_TYPES_H
