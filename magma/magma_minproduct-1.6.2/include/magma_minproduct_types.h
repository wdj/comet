/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015
*/

#ifndef MAGMA_minproduct_TYPES_H
#define MAGMA_minproduct_TYPES_H

#include <stdint.h>
#include <assert.h>


// each implementation of MAGMA_minproduct defines HAVE_* appropriately.
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
// Similar to magma_minproduct_int_t we declare magma_minproduct_index_t used for row/column indices in sparse
#if defined(MAGMA_minproduct_ILP64) || defined(MKL_ILP64)
//typedef int64_t magma_minproduct_int_t;
typedef long long int magma_minproduct_int_t;  // MKL uses long long int, not int64_t
#else
typedef int magma_minproduct_int_t;
#endif

typedef int magma_minproduct_index_t;

// Define new type that the precision generator will not change (matches PLASMA)
typedef double real_Double_t;


// ========================================
// define types specific to implementation (CUDA, OpenCL, MIC)
// define macros to deal with complex numbers
#if defined(HAVE_CUBLAS)
    #ifndef CUBLAS_V2_H_
    #include <cublas.h>
    #endif
    
    typedef cudaStream_t   magma_minproduct_queue_t;
    typedef cudaEvent_t    magma_minproduct_event_t;
    typedef int            magma_minproduct_device_t;
    
    typedef cuDoubleComplex magma_minproductDoubleComplex;
    typedef cuFloatComplex  magma_minproductFloatComplex;
    
    #define MAGMA_minproduct_Z_MAKE(r,i)     make_cuDoubleComplex(r, i)
    #define MAGMA_minproduct_Z_REAL(a)       (a).x
    #define MAGMA_minproduct_Z_IMAG(a)       (a).y
    #define MAGMA_minproduct_Z_ADD(a, b)     cuCadd(a, b)
    #define MAGMA_minproduct_Z_SUB(a, b)     cuCsub(a, b)
    #define MAGMA_minproduct_Z_MUL(a, b)     cuCmul(a, b)
    #define MAGMA_minproduct_Z_DIV(a, b)     cuCdiv(a, b)
    #define MAGMA_minproduct_Z_ABS(a)        cuCabs(a)
    #define MAGMA_minproduct_Z_ABS1(a)       (fabs((a).x) + fabs((a).y))
    #define MAGMA_minproduct_Z_CNJG(a)       cuConj(a)
    
    #define MAGMA_minproduct_C_MAKE(r,i)     make_cuFloatComplex(r, i)
    #define MAGMA_minproduct_C_REAL(a)       (a).x
    #define MAGMA_minproduct_C_IMAG(a)       (a).y
    #define MAGMA_minproduct_C_ADD(a, b)     cuCaddf(a, b)
    #define MAGMA_minproduct_C_SUB(a, b)     cuCsubf(a, b)
    #define MAGMA_minproduct_C_MUL(a, b)     cuCmulf(a, b)
    #define MAGMA_minproduct_C_DIV(a, b)     cuCdivf(a, b)
    #define MAGMA_minproduct_C_ABS(a)        cuCabsf(a)
    #define MAGMA_minproduct_C_ABS1(a)       (fabsf((a).x) + fabsf((a).y))
    #define MAGMA_minproduct_C_CNJG(a)       cuConjf(a)
    
#elif defined(HAVE_clAmdBlas)
    #include <clAmdBlas.h>
    
    typedef cl_command_queue  magma_minproduct_queue_t;
    typedef cl_event          magma_minproduct_event_t;
    typedef cl_device_id      magma_minproduct_device_t;
    
    typedef DoubleComplex magma_minproductDoubleComplex;
    typedef FloatComplex  magma_minproductFloatComplex;
    
    #define MAGMA_minproduct_Z_MAKE(r,i)     doubleComplex(r,i)
    #define MAGMA_minproduct_Z_REAL(a)       (a).x
    #define MAGMA_minproduct_Z_IMAG(a)       (a).y
    #define MAGMA_minproduct_Z_ADD(a, b)     MAGMA_minproduct_Z_MAKE((a).x+(b).x, (a).y+(b).y)
    #define MAGMA_minproduct_Z_SUB(a, b)     MAGMA_minproduct_Z_MAKE((a).x-(b).x, (a).y-(b).y)
    #define MAGMA_minproduct_Z_ABS(a)        magma_minproduct_cabs(a)
    #define MAGMA_minproduct_Z_ABS1(a)       (fabs((a).x) + fabs((a).y))
    #define MAGMA_minproduct_Z_CNJG(a)       MAGMA_minproduct_Z_MAKE((a).x, -(a).y)
    
    #define MAGMA_minproduct_C_MAKE(r,i)     floatComplex(r,i)
    #define MAGMA_minproduct_C_REAL(a)       (a).x
    #define MAGMA_minproduct_C_IMAG(a)       (a).y
    #define MAGMA_minproduct_C_ADD(a, b)     MAGMA_minproduct_C_MAKE((a).x+(b).x, (a).y+(b).y)
    #define MAGMA_minproduct_C_SUB(a, b)     MAGMA_minproduct_C_MAKE((a).x-(b).x, (a).y-(b).y)
    #define MAGMA_minproduct_C_ABS(a)        magma_minproduct_cabsf(a)
    #define MAGMA_minproduct_C_ABS1(a)       (fabsf((a).x) + fabsf((a).y))
    #define MAGMA_minproduct_C_CNJG(a)       MAGMA_minproduct_C_MAKE((a).x, -(a).y)

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

    typedef int   magma_minproduct_queue_t;
    typedef int   magma_minproduct_event_t;
    typedef int   magma_minproduct_device_t;

    #include <complex>
    typedef std::complex<float>   magma_minproductFloatComplex;
    typedef std::complex<double>  magma_minproductDoubleComplex;

    #define MAGMA_minproduct_Z_MAKE(r, i)    std::complex<double>(r,i)
    #define MAGMA_minproduct_Z_REAL(x)       (x).real()
    #define MAGMA_minproduct_Z_IMAG(x)       (x).imag()
    #define MAGMA_minproduct_Z_SET2REAL(a,r) { (a).real() = (r);   (a).imag() = 0.0; }
    #define MAGMA_minproduct_Z_ADD(a, b)     ((a)+(b))
    #define MAGMA_minproduct_Z_SUB(a, b)     ((a)-(b))
    #define MAGMA_minproduct_Z_ABS(a)        abs(a)
    #define MAGMA_minproduct_Z_MUL(a, b)     ((a)*(b))
    #define MAGMA_minproduct_Z_DIV(a, b)     ((a)/(b))
    #define MAGMA_minproduct_Z_ABS1(a)       (fabs((a).real()) + fabs((a).imag())) 
    #define MAGMA_minproduct_Z_CNJG(a)       conj(a)
    #define MAGMA_minproduct_Z_DSCALE(v,t,s) ((v) = (t)/(s))

    #define MAGMA_minproduct_C_MAKE(r, i)    std::complex<float> (r,i)
    #define MAGMA_minproduct_C_REAL(x)       (x).real()
    #define MAGMA_minproduct_C_IMAG(x)       (x).imag()
    #define MAGMA_minproduct_C_SET2REAL(a,r) { (a).real() = (r);   (a).imag() = 0.0; }
    #define MAGMA_minproduct_C_ADD(a, b)     ((a)+(b))
    #define MAGMA_minproduct_C_SUB(a, b)     ((a)-(b))
    #define MAGMA_minproduct_C_ABS(a)        abs(a)
    #define MAGMA_minproduct_C_MUL(a, b)     ((a)*(b))
    #define MAGMA_minproduct_C_DIV(a, b)     ((a)/(b))
    #define MAGMA_minproduct_C_ABS1(a)       (fabs((a).real()) + fabs((a).imag()))
    #define MAGMA_minproduct_C_CNJG(a)       conj(a)
    #define MAGMA_minproduct_C_SSCALE(v,t,s) ((v) = (t)/(s))
#else
    #error "One of HAVE_CUBLAS, HAVE_clAmdBlas, or HAVE_MIC must be defined. For example, add -DHAVE_CUBLAS to CFLAGS, or #define HAVE_CUBLAS before #include <magma_minproduct.h>. In MAGMA_minproduct, this happens in Makefile.internal."
#endif

#define MAGMA_minproduct_Z_EQUAL(a,b)        (MAGMA_minproduct_Z_REAL(a)==MAGMA_minproduct_Z_REAL(b) && MAGMA_minproduct_Z_IMAG(a)==MAGMA_minproduct_Z_IMAG(b))
#define MAGMA_minproduct_Z_NEGATE(a)         MAGMA_minproduct_Z_MAKE( -MAGMA_minproduct_Z_REAL(a), -MAGMA_minproduct_Z_IMAG(a))

#define MAGMA_minproduct_C_EQUAL(a,b)        (MAGMA_minproduct_C_REAL(a)==MAGMA_minproduct_C_REAL(b) && MAGMA_minproduct_C_IMAG(a)==MAGMA_minproduct_C_IMAG(b))
#define MAGMA_minproduct_C_NEGATE(a)         MAGMA_minproduct_C_MAKE( -MAGMA_minproduct_C_REAL(a), -MAGMA_minproduct_C_IMAG(a))

#define MAGMA_minproduct_D_MAKE(r,i)         (r)
#define MAGMA_minproduct_D_REAL(x)           (x)
#define MAGMA_minproduct_D_IMAG(x)           (0.0)
#define MAGMA_minproduct_D_SET2REAL(a,r)     (a) = (r)
#define MAGMA_minproduct_D_ADD(a, b)         ((a) + (b))
#define MAGMA_minproduct_D_SUB(a, b)         ((a) - (b))
#define MAGMA_minproduct_D_MUL(a, b)         ((a) * (b))
#define MAGMA_minproduct_D_DIV(a, b)         ((a) / (b))
#define MAGMA_minproduct_D_ABS(a)            ((a)>0 ? (a) : -(a))
#define MAGMA_minproduct_D_ABS1(a)           ((a)>0 ? (a) : -(a))
#define MAGMA_minproduct_D_CNJG(a)           (a)
#define MAGMA_minproduct_D_EQUAL(a,b)        ((a) == (b))
#define MAGMA_minproduct_D_NEGATE(a)         (-a)
#define MAGMA_minproduct_D_DSCALE(v, t, s)   (v) = (t)/(s)

#define MAGMA_minproduct_S_MAKE(r,i)         (r)
#define MAGMA_minproduct_S_REAL(x)           (x)
#define MAGMA_minproduct_S_IMAG(x)           (0.0)
#define MAGMA_minproduct_S_SET2REAL(a,r)     (a) = (r)
#define MAGMA_minproduct_S_ADD(a, b)         ((a) + (b))
#define MAGMA_minproduct_S_SUB(a, b)         ((a) - (b))
#define MAGMA_minproduct_S_MUL(a, b)         ((a) * (b))
#define MAGMA_minproduct_S_DIV(a, b)         ((a) / (b))
#define MAGMA_minproduct_S_ABS(a)            ((a)>0 ? (a) : -(a))
#define MAGMA_minproduct_S_ABS1(a)           ((a)>0 ? (a) : -(a))
#define MAGMA_minproduct_S_CNJG(a)           (a)
#define MAGMA_minproduct_S_EQUAL(a,b)        ((a) == (b))
#define MAGMA_minproduct_S_NEGATE(a)         (-a)
#define MAGMA_minproduct_S_SSCALE(v, t, s)   (v) = (t)/(s)

#define MAGMA_minproduct_Z_ZERO              MAGMA_minproduct_Z_MAKE( 0.0, 0.0)
#define MAGMA_minproduct_Z_ONE               MAGMA_minproduct_Z_MAKE( 1.0, 0.0)
#define MAGMA_minproduct_Z_HALF              MAGMA_minproduct_Z_MAKE( 0.5, 0.0)
#define MAGMA_minproduct_Z_NEG_ONE           MAGMA_minproduct_Z_MAKE(-1.0, 0.0)
#define MAGMA_minproduct_Z_NEG_HALF          MAGMA_minproduct_Z_MAKE(-0.5, 0.0)

#define MAGMA_minproduct_C_ZERO              MAGMA_minproduct_C_MAKE( 0.0, 0.0)
#define MAGMA_minproduct_C_ONE               MAGMA_minproduct_C_MAKE( 1.0, 0.0)
#define MAGMA_minproduct_C_HALF              MAGMA_minproduct_C_MAKE( 0.5, 0.0)
#define MAGMA_minproduct_C_NEG_ONE           MAGMA_minproduct_C_MAKE(-1.0, 0.0)
#define MAGMA_minproduct_C_NEG_HALF          MAGMA_minproduct_C_MAKE(-0.5, 0.0)

#define MAGMA_minproduct_D_ZERO              ( 0.0)
#define MAGMA_minproduct_D_ONE               ( 1.0)
#define MAGMA_minproduct_D_HALF              ( 0.5)
#define MAGMA_minproduct_D_NEG_ONE           (-1.0)
#define MAGMA_minproduct_D_NEG_HALF          (-0.5)

#define MAGMA_minproduct_S_ZERO              ( 0.0)
#define MAGMA_minproduct_S_ONE               ( 1.0)
#define MAGMA_minproduct_S_HALF              ( 0.5)
#define MAGMA_minproduct_S_NEG_ONE           (-1.0)
#define MAGMA_minproduct_S_NEG_HALF          (-0.5)

#ifndef CBLAS_SADDR
#define CBLAS_SADDR(a)  &(a)
#endif

#if defined(HAVE_clAmdBlas)
    // OpenCL uses opaque memory references on GPU
    typedef cl_mem magma_minproduct_ptr;
    typedef cl_mem magma_minproductInt_ptr;
    typedef cl_mem magma_minproductIndex_ptr;
    typedef cl_mem magma_minproductFloat_ptr;
    typedef cl_mem magma_minproductDouble_ptr;
    typedef cl_mem magma_minproductFloatComplex_ptr;
    typedef cl_mem magma_minproductDoubleComplex_ptr;
    
    typedef cl_mem magma_minproduct_const_ptr;
    typedef cl_mem magma_minproductInt_const_ptr;
    typedef cl_mem magma_minproductIndex_const_ptr;
    typedef cl_mem magma_minproductFloat_const_ptr;
    typedef cl_mem magma_minproductDouble_const_ptr;
    typedef cl_mem magma_minproductFloatComplex_const_ptr;
    typedef cl_mem magma_minproductDoubleComplex_const_ptr;
#else
    // MIC and CUDA use regular pointers on GPU
    typedef void               *magma_minproduct_ptr;
    typedef magma_minproduct_int_t        *magma_minproductInt_ptr;
    typedef magma_minproduct_index_t      *magma_minproductIndex_ptr;
    typedef float              *magma_minproductFloat_ptr;
    typedef double             *magma_minproductDouble_ptr;
    typedef magma_minproductFloatComplex  *magma_minproductFloatComplex_ptr;
    typedef magma_minproductDoubleComplex *magma_minproductDoubleComplex_ptr;
    
    typedef void               const *magma_minproduct_const_ptr;
    typedef magma_minproduct_int_t        const *magma_minproductInt_const_ptr;
    typedef magma_minproduct_index_t      const *magma_minproductIndex_const_ptr;
    typedef float              const *magma_minproductFloat_const_ptr;
    typedef double             const *magma_minproductDouble_const_ptr;
    typedef magma_minproductFloatComplex  const *magma_minproductFloatComplex_const_ptr;
    typedef magma_minproductDoubleComplex const *magma_minproductDoubleComplex_const_ptr;
#endif


// ========================================
// MAGMA_minproduct constants

// ----------------------------------------
#define MAGMA_minproduct_VERSION_MAJOR 1
#define MAGMA_minproduct_VERSION_MINOR 6
#define MAGMA_minproduct_VERSION_MICRO 2

// stage is "svn", "beta#", "rc#" (release candidate), or blank ("") for final release
#define MAGMA_minproduct_VERSION_STAGE ""

#define Magma_minproductMaxGPUs 8
#define Magma_minproductMaxDevices 8

// ----------------------------------------
// Return codes
// LAPACK argument errors are < 0 but > MAGMA_minproduct_ERR.
// MAGMA_minproduct errors are < MAGMA_minproduct_ERR.
#define MAGMA_minproduct_SUCCESS               0
#define MAGMA_minproduct_ERR                  -100
#define MAGMA_minproduct_ERR_NOT_INITIALIZED  -101
#define MAGMA_minproduct_ERR_REINITIALIZED    -102
#define MAGMA_minproduct_ERR_NOT_SUPPORTED    -103
#define MAGMA_minproduct_ERR_ILLEGAL_VALUE    -104
#define MAGMA_minproduct_ERR_NOT_FOUND        -105
#define MAGMA_minproduct_ERR_ALLOCATION       -106
#define MAGMA_minproduct_ERR_INTERNAL_LIMIT   -107
#define MAGMA_minproduct_ERR_UNALLOCATED      -108
#define MAGMA_minproduct_ERR_FILESYSTEM       -109
#define MAGMA_minproduct_ERR_UNEXPECTED       -110
#define MAGMA_minproduct_ERR_SEQUENCE_FLUSHED -111
#define MAGMA_minproduct_ERR_HOST_ALLOC       -112
#define MAGMA_minproduct_ERR_DEVICE_ALLOC     -113
#define MAGMA_minproduct_ERR_CUDASTREAM       -114
#define MAGMA_minproduct_ERR_INVALID_PTR      -115
#define MAGMA_minproduct_ERR_UNKNOWN          -116
#define MAGMA_minproduct_ERR_NOT_IMPLEMENTED  -117

// some sparse-iter errors
#define MAGMA_minproduct_SLOW_CONVERGENCE     -201
#define MAGMA_minproduct_DIVERGENCE           -202
#define MAGMA_minproduct_NONSPD               -203
#define MAGMA_minproduct_ERR_BADPRECOND       -204

// When adding error codes, please add to interface_cuda/error.cpp

// map cusparse errors to magma_minproduct errors
#define MAGMA_minproduct_ERR_CUSPARSE_NOT_INITIALIZED            -3001
#define MAGMA_minproduct_ERR_CUSPARSE_ALLOC_FAILED               -3002
#define MAGMA_minproduct_ERR_CUSPARSE_INVALID_VALUE              -3003
#define MAGMA_minproduct_ERR_CUSPARSE_ARCH_MISMATCH              -3004
#define MAGMA_minproduct_ERR_CUSPARSE_MAPPING_ERROR              -3005
#define MAGMA_minproduct_ERR_CUSPARSE_EXECUTION_FAILED           -3006
#define MAGMA_minproduct_ERR_CUSPARSE_INTERNAL_ERROR             -3007
#define MAGMA_minproduct_ERR_CUSPARSE_MATRIX_TYPE_NOT_SUPPORTED  -3008
#define MAGMA_minproduct_ERR_CUSPARSE_ZERO_PIVOT                 -3009


// ----------------------------------------
// parameter constants
// numbering is consistent with CBLAS and PLASMA; see plasma/include/plasma.h
// also with lapack_cwrapper/include/lapack_enum.h
typedef enum {
    Magma_minproductFalse         = 0,
    Magma_minproductTrue          = 1
} magma_minproduct_bool_t;

typedef enum {
    Magma_minproductRowMajor      = 101,
    Magma_minproductColMajor      = 102
} magma_minproduct_order_t;

// Magma_minproduct_ConjTrans is an alias for those rare occasions (zlarfb, zun*, zher*k)
// where we want Magma_minproduct_ConjTrans to convert to Magma_minproductTrans in precision generation.
typedef enum {
    Magma_minproductNoTrans       = 111,
    Magma_minproductTrans         = 112,
    Magma_minproductConjTrans     = 113,
    Magma_minproduct_ConjTrans    = Magma_minproductConjTrans
} magma_minproduct_trans_t;

typedef enum {
    Magma_minproductUpper         = 121,
    Magma_minproductLower         = 122,
    Magma_minproductUpperLower    = 123,
    Magma_minproductFull          = 123,  /* lascl, laset */
    Magma_minproductHessenberg    = 124   /* lascl */
} magma_minproduct_uplo_t;

typedef magma_minproduct_uplo_t magma_minproduct_type_t;  /* lascl */

typedef enum {
    Magma_minproductNonUnit       = 131,
    Magma_minproductUnit          = 132
} magma_minproduct_diag_t;

typedef enum {
    Magma_minproductLeft          = 141,
    Magma_minproductRight         = 142,
    Magma_minproductBothSides     = 143   /* trevc */
} magma_minproduct_side_t;

typedef enum {
    Magma_minproductOneNorm       = 171,  /* lange, lanhe */
    Magma_minproductRealOneNorm   = 172,
    Magma_minproductTwoNorm       = 173,
    Magma_minproductFrobeniusNorm = 174,
    Magma_minproductInfNorm       = 175,
    Magma_minproductRealInfNorm   = 176,
    Magma_minproductMaxNorm       = 177,
    Magma_minproductRealMaxNorm   = 178
} magma_minproduct_norm_t;

typedef enum {
    Magma_minproductDistUniform   = 201,  /* latms */
    Magma_minproductDistSymmetric = 202,
    Magma_minproductDistNormal    = 203
} magma_minproduct_dist_t;

typedef enum {
    Magma_minproductHermGeev      = 241,  /* latms */
    Magma_minproductHermPoev      = 242,
    Magma_minproductNonsymPosv    = 243,
    Magma_minproductSymPosv       = 244
} magma_minproduct_sym_t;

typedef enum {
    Magma_minproductNoPacking     = 291,  /* latms */
    Magma_minproductPackSubdiag   = 292,
    Magma_minproductPackSupdiag   = 293,
    Magma_minproductPackColumn    = 294,
    Magma_minproductPackRow       = 295,
    Magma_minproductPackLowerBand = 296,
    Magma_minproductPackUpeprBand = 297,
    Magma_minproductPackAll       = 298
} magma_minproduct_pack_t;

typedef enum {
    Magma_minproductNoVec         = 301,  /* geev, syev, gesvd */
    Magma_minproductVec           = 302,  /* geev, syev */
    Magma_minproductIVec          = 303,  /* stedc */
    Magma_minproductAllVec        = 304,  /* gesvd, trevc */
    Magma_minproductSomeVec       = 305,  /* gesvd, trevc */
    Magma_minproductOverwriteVec  = 306,  /* gesvd */
    Magma_minproductBacktransVec  = 307   /* trevc */
} magma_minproduct_vec_t;

typedef enum {
    Magma_minproductRangeAll      = 311,  /* syevx, etc. */
    Magma_minproductRangeV        = 312,
    Magma_minproductRangeI        = 313
} magma_minproduct_range_t;

typedef enum {
    Magma_minproductQ             = 322,  /* unmbr, ungbr */
    Magma_minproductP             = 323
} magma_minproduct_vect_t;

typedef enum {
    Magma_minproductForward       = 391,  /* larfb */
    Magma_minproductBackward      = 392
} magma_minproduct_direct_t;

typedef enum {
    Magma_minproductColumnwise    = 401,  /* larfb */
    Magma_minproductRowwise       = 402
} magma_minproduct_storev_t;

// --------------------
// sparse
typedef enum {
    Magma_minproduct_CSR          = 411,
    Magma_minproduct_ELLPACKT     = 412,
    Magma_minproduct_ELL          = 413,
    Magma_minproduct_DENSE        = 414,
    Magma_minproduct_BCSR         = 415,
    Magma_minproduct_CSC          = 416,
    Magma_minproduct_HYB          = 417,
    Magma_minproduct_COO          = 418,
    Magma_minproduct_ELLRT        = 419,
    Magma_minproduct_SPMVFUNCTION = 420,
    Magma_minproduct_SELLP        = 421,
    Magma_minproduct_ELLD         = 422,

    Magma_minproduct_CSRD         = 424,
    Magma_minproduct_CSRL         = 427,
    Magma_minproduct_CSRU         = 428,
    Magma_minproduct_CSRCOO       = 429
} magma_minproduct_storage_t;


typedef enum {
    Magma_minproduct_CG           = 431,
    Magma_minproduct_CGMERGE      = 432,
    Magma_minproduct_GMRES        = 433,
    Magma_minproduct_BICGSTAB     = 434,
  Magma_minproduct_BICGSTABMERGE  = 435,
  Magma_minproduct_BICGSTABMERGE2 = 436,
    Magma_minproduct_JACOBI       = 437,
    Magma_minproduct_GS           = 438,
    Magma_minproduct_ITERREF      = 439,
    Magma_minproduct_BCSRLU       = 440,
    Magma_minproduct_PCG          = 441,
    Magma_minproduct_PGMRES       = 442,
    Magma_minproduct_PBICGSTAB    = 443,
    Magma_minproduct_PASTIX       = 444,
    Magma_minproduct_ILU          = 445,
    Magma_minproduct_ICC          = 446,
    Magma_minproduct_AILU         = 447,
    Magma_minproduct_AICC         = 448,
    Magma_minproduct_BAITER       = 449,
    Magma_minproduct_LOBPCG       = 450,
    Magma_minproduct_NONE         = 451,
    Magma_minproduct_FUNCTION     = 452
} magma_minproduct_solver_type;

typedef enum {
    Magma_minproduct_CGS          = 461,
    Magma_minproduct_FUSED_CGS    = 462,
    Magma_minproduct_MGS          = 463
} magma_minproduct_ortho_t;

typedef enum {
    Magma_minproduct_CPU          = 471,
    Magma_minproduct_DEV          = 472
} magma_minproduct_location_t;

typedef enum {
    Magma_minproduct_GENERAL      = 481,
    Magma_minproduct_SYMMETRIC    = 482
} magma_minproduct_symmetry_t;

typedef enum {
    Magma_minproduct_ORDERED      = 491,
    Magma_minproduct_DIAGFIRST    = 492,
    Magma_minproduct_UNITY        = 493,
    Magma_minproduct_VALUE        = 494
} magma_minproduct_diagorder_t;

typedef enum {
    Magma_minproduct_DCOMPLEX     = 501,
    Magma_minproduct_FCOMPLEX     = 502,
    Magma_minproduct_DOUBLE       = 503,
    Magma_minproduct_FLOAT        = 504
} magma_minproduct_precision;

typedef enum {
    Magma_minproduct_NOSCALE      = 511,
    Magma_minproduct_UNITROW      = 512,
    Magma_minproduct_UNITDIAG     = 513
} magma_minproduct_scale_t;

typedef enum {
    Magma_minproduct_FULL         = 521,
    Magma_minproduct_LOWER        = 522,
    Magma_minproduct_UPPER        = 523
} magma_minproduct_fillmode_t;



// When adding constants, remember to do these steps as appropriate:
// 1)  add magma_minproduct_xxxx_const()  converter below and in control/constants.cpp
// 2a) add to magma_minproduct2lapack_constants[] in control/constants.cpp
// 2b) update min & max here, which are used to check bounds for magma_minproduct2lapack_constants[]
// 2c) add lapack_xxxx_const() converter below and in control/constants.cpp
#define Magma_minproduct2lapack_Min  Magma_minproductFalse     // 0
#define Magma_minproduct2lapack_Max  Magma_minproductRowwise   // 402


// ----------------------------------------
// string constants for calling Fortran BLAS and LAPACK
// todo: use translators instead? lapack_const( Magma_minproductUpper )
#define Magma_minproductRowMajorStr      "Row"
#define Magma_minproductColMajorStr      "Col"

#define Magma_minproductNoTransStr       "NoTrans"
#define Magma_minproductTransStr         "Trans"
#define Magma_minproductConjTransStr     "ConjTrans"
#define Magma_minproduct_ConjTransStr    "ConjTrans"

#define Magma_minproductUpperStr         "Upper"
#define Magma_minproductLowerStr         "Lower"
#define Magma_minproductUpperLowerStr    "Full"
#define Magma_minproductFullStr          "Full"

#define Magma_minproductNonUnitStr       "NonUnit"
#define Magma_minproductUnitStr          "Unit"

#define Magma_minproductLeftStr          "Left"
#define Magma_minproductRightStr         "Right"
#define Magma_minproductBothSidesStr     "Both"

#define Magma_minproductOneNormStr       "1"
#define Magma_minproductTwoNormStr       "2"
#define Magma_minproductFrobeniusNormStr "Fro"
#define Magma_minproductInfNormStr       "Inf"
#define Magma_minproductMaxNormStr       "Max"

#define Magma_minproductForwardStr       "Forward"
#define Magma_minproductBackwardStr      "Backward"

#define Magma_minproductColumnwiseStr    "Columnwise"
#define Magma_minproductRowwiseStr       "Rowwise"

#define Magma_minproductNoVecStr         "NoVec"
#define Magma_minproductVecStr           "Vec"
#define Magma_minproductIVecStr          "IVec"
#define Magma_minproductAllVecStr        "All"
#define Magma_minproductSomeVecStr       "Some"
#define Magma_minproductOverwriteVecStr  "Overwrite"


#ifdef __cplusplus
extern "C" {
#endif

// --------------------
// Convert LAPACK character constants to MAGMA_minproduct constants.
// This is a one-to-many mapping, requiring multiple translators
// (e.g., "N" can be NoTrans or NonUnit or NoVec).
magma_minproduct_bool_t   magma_minproduct_bool_const  ( char lapack_char );
magma_minproduct_order_t  magma_minproduct_order_const ( char lapack_char );
magma_minproduct_trans_t  magma_minproduct_trans_const ( char lapack_char );
magma_minproduct_uplo_t   magma_minproduct_uplo_const  ( char lapack_char );
magma_minproduct_diag_t   magma_minproduct_diag_const  ( char lapack_char );
magma_minproduct_side_t   magma_minproduct_side_const  ( char lapack_char );
magma_minproduct_norm_t   magma_minproduct_norm_const  ( char lapack_char );
magma_minproduct_dist_t   magma_minproduct_dist_const  ( char lapack_char );
magma_minproduct_sym_t    magma_minproduct_sym_const   ( char lapack_char );
magma_minproduct_pack_t   magma_minproduct_pack_const  ( char lapack_char );
magma_minproduct_vec_t    magma_minproduct_vec_const   ( char lapack_char );
magma_minproduct_range_t  magma_minproduct_range_const ( char lapack_char );
magma_minproduct_vect_t   magma_minproduct_vect_const  ( char lapack_char );
magma_minproduct_direct_t magma_minproduct_direct_const( char lapack_char );
magma_minproduct_storev_t magma_minproduct_storev_const( char lapack_char );


// --------------------
// Convert MAGMA_minproduct constants to LAPACK(E) constants.
// The generic lapack_const works for all cases, but the specific routines
// (e.g., lapack_trans_const) do better error checking.
const char* lapack_const       ( int            magma_minproduct_const );
const char* lapack_bool_const  ( magma_minproduct_bool_t   magma_minproduct_const );
const char* lapack_order_const ( magma_minproduct_order_t  magma_minproduct_const );
const char* lapack_trans_const ( magma_minproduct_trans_t  magma_minproduct_const );
const char* lapack_uplo_const  ( magma_minproduct_uplo_t   magma_minproduct_const );
const char* lapack_diag_const  ( magma_minproduct_diag_t   magma_minproduct_const );
const char* lapack_side_const  ( magma_minproduct_side_t   magma_minproduct_const );
const char* lapack_norm_const  ( magma_minproduct_norm_t   magma_minproduct_const );
const char* lapack_dist_const  ( magma_minproduct_dist_t   magma_minproduct_const );
const char* lapack_sym_const   ( magma_minproduct_sym_t    magma_minproduct_const );
const char* lapack_pack_const  ( magma_minproduct_pack_t   magma_minproduct_const );
const char* lapack_vec_const   ( magma_minproduct_vec_t    magma_minproduct_const );
const char* lapack_range_const ( magma_minproduct_range_t  magma_minproduct_const );
const char* lapack_vect_const  ( magma_minproduct_vect_t   magma_minproduct_const );
const char* lapack_direct_const( magma_minproduct_direct_t magma_minproduct_const );
const char* lapack_storev_const( magma_minproduct_storev_t magma_minproduct_const );

static inline char lapacke_const       ( int magma_minproduct_const            ) { return *lapack_const       ( magma_minproduct_const ); }
static inline char lapacke_bool_const  ( magma_minproduct_bool_t   magma_minproduct_const ) { return *lapack_bool_const  ( magma_minproduct_const ); }
static inline char lapacke_order_const ( magma_minproduct_order_t  magma_minproduct_const ) { return *lapack_order_const ( magma_minproduct_const ); }
static inline char lapacke_trans_const ( magma_minproduct_trans_t  magma_minproduct_const ) { return *lapack_trans_const ( magma_minproduct_const ); }
static inline char lapacke_uplo_const  ( magma_minproduct_uplo_t   magma_minproduct_const ) { return *lapack_uplo_const  ( magma_minproduct_const ); }
static inline char lapacke_diag_const  ( magma_minproduct_diag_t   magma_minproduct_const ) { return *lapack_diag_const  ( magma_minproduct_const ); }
static inline char lapacke_side_const  ( magma_minproduct_side_t   magma_minproduct_const ) { return *lapack_side_const  ( magma_minproduct_const ); }
static inline char lapacke_norm_const  ( magma_minproduct_norm_t   magma_minproduct_const ) { return *lapack_norm_const  ( magma_minproduct_const ); }
static inline char lapacke_dist_const  ( magma_minproduct_dist_t   magma_minproduct_const ) { return *lapack_dist_const  ( magma_minproduct_const ); }
static inline char lapacke_sym_const   ( magma_minproduct_sym_t    magma_minproduct_const ) { return *lapack_sym_const   ( magma_minproduct_const ); }
static inline char lapacke_pack_const  ( magma_minproduct_pack_t   magma_minproduct_const ) { return *lapack_pack_const  ( magma_minproduct_const ); }
static inline char lapacke_vec_const   ( magma_minproduct_vec_t    magma_minproduct_const ) { return *lapack_vec_const   ( magma_minproduct_const ); }
static inline char lapacke_range_const ( magma_minproduct_range_t  magma_minproduct_const ) { return *lapack_range_const ( magma_minproduct_const ); }
static inline char lapacke_vect_const  ( magma_minproduct_vect_t   magma_minproduct_const ) { return *lapack_vect_const  ( magma_minproduct_const ); }
static inline char lapacke_direct_const( magma_minproduct_direct_t magma_minproduct_const ) { return *lapack_direct_const( magma_minproduct_const ); }
static inline char lapacke_storev_const( magma_minproduct_storev_t magma_minproduct_const ) { return *lapack_storev_const( magma_minproduct_const ); }


// --------------------
// Convert MAGMA_minproduct constants to clAmdBlas constants.
#if defined(HAVE_clAmdBlas)
clAmdBlasOrder       amdblas_order_const( magma_minproduct_order_t order );
clAmdBlasTranspose   amdblas_trans_const( magma_minproduct_trans_t trans );
clAmdBlasUplo        amdblas_uplo_const ( magma_minproduct_uplo_t  uplo  );
clAmdBlasDiag        amdblas_diag_const ( magma_minproduct_diag_t  diag  );
clAmdBlasSide        amdblas_side_const ( magma_minproduct_side_t  side  );
#endif


// --------------------
// Convert MAGMA_minproduct constants to CUBLAS constants.
#if defined(CUBLAS_V2_H_)
cublasOperation_t    cublas_trans_const ( magma_minproduct_trans_t trans );
cublasFillMode_t     cublas_uplo_const  ( magma_minproduct_uplo_t  uplo  );
cublasDiagType_t     cublas_diag_const  ( magma_minproduct_diag_t  diag  );
cublasSideMode_t     cublas_side_const  ( magma_minproduct_side_t  side  );
#endif


// --------------------
// Convert MAGMA_minproduct constants to CBLAS constants.
#if defined(HAVE_CBLAS)
#include <cblas.h>
enum CBLAS_ORDER     cblas_order_const  ( magma_minproduct_order_t order );
enum CBLAS_TRANSPOSE cblas_trans_const  ( magma_minproduct_trans_t trans );
enum CBLAS_UPLO      cblas_uplo_const   ( magma_minproduct_uplo_t  uplo  );
enum CBLAS_DIAG      cblas_diag_const   ( magma_minproduct_diag_t  diag  );
enum CBLAS_SIDE      cblas_side_const   ( magma_minproduct_side_t  side  );
#endif


#ifdef __cplusplus
}
#endif

#endif        //  #ifndef MAGMA_minproduct_TYPES_H
