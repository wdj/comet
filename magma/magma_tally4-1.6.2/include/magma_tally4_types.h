/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015
*/

#ifndef MAGMA_tally4_TYPES_H
#define MAGMA_tally4_TYPES_H

#include <stdint.h>
#include <assert.h>


// each implementation of MAGMA_tally4 defines HAVE_* appropriately.
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
// Similar to magma_tally4_int_t we declare magma_tally4_index_t used for row/column indices in sparse
#if defined(MAGMA_tally4_ILP64) || defined(MKL_ILP64)
//typedef int64_t magma_tally4_int_t;
typedef long long int magma_tally4_int_t;  // MKL uses long long int, not int64_t
#else
typedef int magma_tally4_int_t;
#endif

typedef int magma_tally4_index_t;

// Define new type that the precision generator will not change (matches PLASMA)
typedef double real_Double_t;


// ========================================
// define types specific to implementation (CUDA, OpenCL, MIC)
// define macros to deal with complex numbers
#if defined(HAVE_CUBLAS)
    #ifndef CUBLAS_V2_H_
    #include <cublas.h>
    #endif
    
    typedef cudaStream_t   magma_tally4_queue_t;
    typedef cudaEvent_t    magma_tally4_event_t;
    typedef int            magma_tally4_device_t;
    
    typedef cuDoubleComplex magma_tally4DoubleComplex;
    typedef cuFloatComplex  magma_tally4FloatComplex;
    
    #define MAGMA_tally4_Z_MAKE(r,i)     make_cuDoubleComplex(r, i)
    #define MAGMA_tally4_Z_REAL(a)       (a).x
    #define MAGMA_tally4_Z_IMAG(a)       (a).y
    #define MAGMA_tally4_Z_ADD(a, b)     cuCadd(a, b)
    #define MAGMA_tally4_Z_SUB(a, b)     cuCsub(a, b)
    #define MAGMA_tally4_Z_MUL(a, b)     cuCmul(a, b)
    #define MAGMA_tally4_Z_DIV(a, b)     cuCdiv(a, b)
    #define MAGMA_tally4_Z_ABS(a)        cuCabs(a)
    #define MAGMA_tally4_Z_ABS1(a)       (fabs((a).x) + fabs((a).y))
    #define MAGMA_tally4_Z_CNJG(a)       cuConj(a)
    
    #define MAGMA_tally4_C_MAKE(r,i)     make_cuFloatComplex(r, i)
    #define MAGMA_tally4_C_REAL(a)       (a).x
    #define MAGMA_tally4_C_IMAG(a)       (a).y
    #define MAGMA_tally4_C_ADD(a, b)     cuCaddf(a, b)
    #define MAGMA_tally4_C_SUB(a, b)     cuCsubf(a, b)
    #define MAGMA_tally4_C_MUL(a, b)     cuCmulf(a, b)
    #define MAGMA_tally4_C_DIV(a, b)     cuCdivf(a, b)
    #define MAGMA_tally4_C_ABS(a)        cuCabsf(a)
    #define MAGMA_tally4_C_ABS1(a)       (fabsf((a).x) + fabsf((a).y))
    #define MAGMA_tally4_C_CNJG(a)       cuConjf(a)
    
#elif defined(HAVE_clAmdBlas)
    #include <clAmdBlas.h>
    
    typedef cl_command_queue  magma_tally4_queue_t;
    typedef cl_event          magma_tally4_event_t;
    typedef cl_device_id      magma_tally4_device_t;
    
    typedef DoubleComplex magma_tally4DoubleComplex;
    typedef FloatComplex  magma_tally4FloatComplex;
    
    #define MAGMA_tally4_Z_MAKE(r,i)     doubleComplex(r,i)
    #define MAGMA_tally4_Z_REAL(a)       (a).x
    #define MAGMA_tally4_Z_IMAG(a)       (a).y
    #define MAGMA_tally4_Z_ADD(a, b)     MAGMA_tally4_Z_MAKE((a).x+(b).x, (a).y+(b).y)
    #define MAGMA_tally4_Z_SUB(a, b)     MAGMA_tally4_Z_MAKE((a).x-(b).x, (a).y-(b).y)
    #define MAGMA_tally4_Z_ABS(a)        magma_tally4_cabs(a)
    #define MAGMA_tally4_Z_ABS1(a)       (fabs((a).x) + fabs((a).y))
    #define MAGMA_tally4_Z_CNJG(a)       MAGMA_tally4_Z_MAKE((a).x, -(a).y)
    
    #define MAGMA_tally4_C_MAKE(r,i)     floatComplex(r,i)
    #define MAGMA_tally4_C_REAL(a)       (a).x
    #define MAGMA_tally4_C_IMAG(a)       (a).y
    #define MAGMA_tally4_C_ADD(a, b)     MAGMA_tally4_C_MAKE((a).x+(b).x, (a).y+(b).y)
    #define MAGMA_tally4_C_SUB(a, b)     MAGMA_tally4_C_MAKE((a).x-(b).x, (a).y-(b).y)
    #define MAGMA_tally4_C_ABS(a)        magma_tally4_cabsf(a)
    #define MAGMA_tally4_C_ABS1(a)       (fabsf((a).x) + fabsf((a).y))
    #define MAGMA_tally4_C_CNJG(a)       MAGMA_tally4_C_MAKE((a).x, -(a).y)

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

    typedef int   magma_tally4_queue_t;
    typedef int   magma_tally4_event_t;
    typedef int   magma_tally4_device_t;

    #include <complex>
    typedef std::complex<float>   magma_tally4FloatComplex;
    typedef std::complex<double>  magma_tally4DoubleComplex;

    #define MAGMA_tally4_Z_MAKE(r, i)    std::complex<double>(r,i)
    #define MAGMA_tally4_Z_REAL(x)       (x).real()
    #define MAGMA_tally4_Z_IMAG(x)       (x).imag()
    #define MAGMA_tally4_Z_SET2REAL(a,r) { (a).real() = (r);   (a).imag() = 0.0; }
    #define MAGMA_tally4_Z_ADD(a, b)     ((a)+(b))
    #define MAGMA_tally4_Z_SUB(a, b)     ((a)-(b))
    #define MAGMA_tally4_Z_ABS(a)        abs(a)
    #define MAGMA_tally4_Z_MUL(a, b)     ((a)*(b))
    #define MAGMA_tally4_Z_DIV(a, b)     ((a)/(b))
    #define MAGMA_tally4_Z_ABS1(a)       (fabs((a).real()) + fabs((a).imag())) 
    #define MAGMA_tally4_Z_CNJG(a)       conj(a)
    #define MAGMA_tally4_Z_DSCALE(v,t,s) ((v) = (t)/(s))

    #define MAGMA_tally4_C_MAKE(r, i)    std::complex<float> (r,i)
    #define MAGMA_tally4_C_REAL(x)       (x).real()
    #define MAGMA_tally4_C_IMAG(x)       (x).imag()
    #define MAGMA_tally4_C_SET2REAL(a,r) { (a).real() = (r);   (a).imag() = 0.0; }
    #define MAGMA_tally4_C_ADD(a, b)     ((a)+(b))
    #define MAGMA_tally4_C_SUB(a, b)     ((a)-(b))
    #define MAGMA_tally4_C_ABS(a)        abs(a)
    #define MAGMA_tally4_C_MUL(a, b)     ((a)*(b))
    #define MAGMA_tally4_C_DIV(a, b)     ((a)/(b))
    #define MAGMA_tally4_C_ABS1(a)       (fabs((a).real()) + fabs((a).imag()))
    #define MAGMA_tally4_C_CNJG(a)       conj(a)
    #define MAGMA_tally4_C_SSCALE(v,t,s) ((v) = (t)/(s))
#else
    #error "One of HAVE_CUBLAS, HAVE_clAmdBlas, or HAVE_MIC must be defined. For example, add -DHAVE_CUBLAS to CFLAGS, or #define HAVE_CUBLAS before #include <magma_tally4.h>. In MAGMA_tally4, this happens in Makefile.internal."
#endif

#define MAGMA_tally4_Z_EQUAL(a,b)        (MAGMA_tally4_Z_REAL(a)==MAGMA_tally4_Z_REAL(b) && MAGMA_tally4_Z_IMAG(a)==MAGMA_tally4_Z_IMAG(b))
#define MAGMA_tally4_Z_NEGATE(a)         MAGMA_tally4_Z_MAKE( -MAGMA_tally4_Z_REAL(a), -MAGMA_tally4_Z_IMAG(a))

#define MAGMA_tally4_C_EQUAL(a,b)        (MAGMA_tally4_C_REAL(a)==MAGMA_tally4_C_REAL(b) && MAGMA_tally4_C_IMAG(a)==MAGMA_tally4_C_IMAG(b))
#define MAGMA_tally4_C_NEGATE(a)         MAGMA_tally4_C_MAKE( -MAGMA_tally4_C_REAL(a), -MAGMA_tally4_C_IMAG(a))

#define MAGMA_tally4_D_MAKE(r,i)         (r)
#define MAGMA_tally4_D_REAL(x)           (x)
#define MAGMA_tally4_D_IMAG(x)           (0.0)
#define MAGMA_tally4_D_SET2REAL(a,r)     (a) = (r)
#define MAGMA_tally4_D_ADD(a, b)         ((a) + (b))
#define MAGMA_tally4_D_SUB(a, b)         ((a) - (b))
#define MAGMA_tally4_D_MUL(a, b)         ((a) * (b))
#define MAGMA_tally4_D_DIV(a, b)         ((a) / (b))
#define MAGMA_tally4_D_ABS(a)            ((a)>0 ? (a) : -(a))
#define MAGMA_tally4_D_ABS1(a)           ((a)>0 ? (a) : -(a))
#define MAGMA_tally4_D_CNJG(a)           (a)
#define MAGMA_tally4_D_EQUAL(a,b)        ((a) == (b))
#define MAGMA_tally4_D_NEGATE(a)         (-a)
#define MAGMA_tally4_D_DSCALE(v, t, s)   (v) = (t)/(s)

#define MAGMA_tally4_S_MAKE(r,i)         (r)
#define MAGMA_tally4_S_REAL(x)           (x)
#define MAGMA_tally4_S_IMAG(x)           (0.0)
#define MAGMA_tally4_S_SET2REAL(a,r)     (a) = (r)
#define MAGMA_tally4_S_ADD(a, b)         ((a) + (b))
#define MAGMA_tally4_S_SUB(a, b)         ((a) - (b))
#define MAGMA_tally4_S_MUL(a, b)         ((a) * (b))
#define MAGMA_tally4_S_DIV(a, b)         ((a) / (b))
#define MAGMA_tally4_S_ABS(a)            ((a)>0 ? (a) : -(a))
#define MAGMA_tally4_S_ABS1(a)           ((a)>0 ? (a) : -(a))
#define MAGMA_tally4_S_CNJG(a)           (a)
#define MAGMA_tally4_S_EQUAL(a,b)        ((a) == (b))
#define MAGMA_tally4_S_NEGATE(a)         (-a)
#define MAGMA_tally4_S_SSCALE(v, t, s)   (v) = (t)/(s)

#define MAGMA_tally4_Z_ZERO              MAGMA_tally4_Z_MAKE( 0.0, 0.0)
#define MAGMA_tally4_Z_ONE               MAGMA_tally4_Z_MAKE( 1.0, 0.0)
#define MAGMA_tally4_Z_HALF              MAGMA_tally4_Z_MAKE( 0.5, 0.0)
#define MAGMA_tally4_Z_NEG_ONE           MAGMA_tally4_Z_MAKE(-1.0, 0.0)
#define MAGMA_tally4_Z_NEG_HALF          MAGMA_tally4_Z_MAKE(-0.5, 0.0)

#define MAGMA_tally4_C_ZERO              MAGMA_tally4_C_MAKE( 0.0, 0.0)
#define MAGMA_tally4_C_ONE               MAGMA_tally4_C_MAKE( 1.0, 0.0)
#define MAGMA_tally4_C_HALF              MAGMA_tally4_C_MAKE( 0.5, 0.0)
#define MAGMA_tally4_C_NEG_ONE           MAGMA_tally4_C_MAKE(-1.0, 0.0)
#define MAGMA_tally4_C_NEG_HALF          MAGMA_tally4_C_MAKE(-0.5, 0.0)

#define MAGMA_tally4_D_ZERO              ( 0.0)
#define MAGMA_tally4_D_ONE               ( 1.0)
#define MAGMA_tally4_D_HALF              ( 0.5)
#define MAGMA_tally4_D_NEG_ONE           (-1.0)
#define MAGMA_tally4_D_NEG_HALF          (-0.5)

#define MAGMA_tally4_S_ZERO              ( 0.0)
#define MAGMA_tally4_S_ONE               ( 1.0)
#define MAGMA_tally4_S_HALF              ( 0.5)
#define MAGMA_tally4_S_NEG_ONE           (-1.0)
#define MAGMA_tally4_S_NEG_HALF          (-0.5)

#ifndef CBLAS_SADDR
#define CBLAS_SADDR(a)  &(a)
#endif

#if defined(HAVE_clAmdBlas)
    // OpenCL uses opaque memory references on GPU
    typedef cl_mem magma_tally4_ptr;
    typedef cl_mem magma_tally4Int_ptr;
    typedef cl_mem magma_tally4Index_ptr;
    typedef cl_mem magma_tally4Float_ptr;
    typedef cl_mem magma_tally4Double_ptr;
    typedef cl_mem magma_tally4FloatComplex_ptr;
    typedef cl_mem magma_tally4DoubleComplex_ptr;
    
    typedef cl_mem magma_tally4_const_ptr;
    typedef cl_mem magma_tally4Int_const_ptr;
    typedef cl_mem magma_tally4Index_const_ptr;
    typedef cl_mem magma_tally4Float_const_ptr;
    typedef cl_mem magma_tally4Double_const_ptr;
    typedef cl_mem magma_tally4FloatComplex_const_ptr;
    typedef cl_mem magma_tally4DoubleComplex_const_ptr;
#else
    // MIC and CUDA use regular pointers on GPU
    typedef void               *magma_tally4_ptr;
    typedef magma_tally4_int_t        *magma_tally4Int_ptr;
    typedef magma_tally4_index_t      *magma_tally4Index_ptr;
    typedef float              *magma_tally4Float_ptr;
    typedef double             *magma_tally4Double_ptr;
    typedef magma_tally4FloatComplex  *magma_tally4FloatComplex_ptr;
    typedef magma_tally4DoubleComplex *magma_tally4DoubleComplex_ptr;
    
    typedef void               const *magma_tally4_const_ptr;
    typedef magma_tally4_int_t        const *magma_tally4Int_const_ptr;
    typedef magma_tally4_index_t      const *magma_tally4Index_const_ptr;
    typedef float              const *magma_tally4Float_const_ptr;
    typedef double             const *magma_tally4Double_const_ptr;
    typedef magma_tally4FloatComplex  const *magma_tally4FloatComplex_const_ptr;
    typedef magma_tally4DoubleComplex const *magma_tally4DoubleComplex_const_ptr;
#endif


// ========================================
// MAGMA_tally4 constants

// ----------------------------------------
#define MAGMA_tally4_VERSION_MAJOR 1
#define MAGMA_tally4_VERSION_MINOR 6
#define MAGMA_tally4_VERSION_MICRO 2

// stage is "svn", "beta#", "rc#" (release candidate), or blank ("") for final release
#define MAGMA_tally4_VERSION_STAGE ""

#define Magma_tally4MaxGPUs 8
#define Magma_tally4MaxDevices 8

// ----------------------------------------
// Return codes
// LAPACK argument errors are < 0 but > MAGMA_tally4_ERR.
// MAGMA_tally4 errors are < MAGMA_tally4_ERR.
#define MAGMA_tally4_SUCCESS               0
#define MAGMA_tally4_ERR                  -100
#define MAGMA_tally4_ERR_NOT_INITIALIZED  -101
#define MAGMA_tally4_ERR_REINITIALIZED    -102
#define MAGMA_tally4_ERR_NOT_SUPPORTED    -103
#define MAGMA_tally4_ERR_ILLEGAL_VALUE    -104
#define MAGMA_tally4_ERR_NOT_FOUND        -105
#define MAGMA_tally4_ERR_ALLOCATION       -106
#define MAGMA_tally4_ERR_INTERNAL_LIMIT   -107
#define MAGMA_tally4_ERR_UNALLOCATED      -108
#define MAGMA_tally4_ERR_FILESYSTEM       -109
#define MAGMA_tally4_ERR_UNEXPECTED       -110
#define MAGMA_tally4_ERR_SEQUENCE_FLUSHED -111
#define MAGMA_tally4_ERR_HOST_ALLOC       -112
#define MAGMA_tally4_ERR_DEVICE_ALLOC     -113
#define MAGMA_tally4_ERR_CUDASTREAM       -114
#define MAGMA_tally4_ERR_INVALID_PTR      -115
#define MAGMA_tally4_ERR_UNKNOWN          -116
#define MAGMA_tally4_ERR_NOT_IMPLEMENTED  -117

// some sparse-iter errors
#define MAGMA_tally4_SLOW_CONVERGENCE     -201
#define MAGMA_tally4_DIVERGENCE           -202
#define MAGMA_tally4_NONSPD               -203
#define MAGMA_tally4_ERR_BADPRECOND       -204

// When adding error codes, please add to interface_cuda/error.cpp

// map cusparse errors to magma_tally4 errors
#define MAGMA_tally4_ERR_CUSPARSE_NOT_INITIALIZED            -3001
#define MAGMA_tally4_ERR_CUSPARSE_ALLOC_FAILED               -3002
#define MAGMA_tally4_ERR_CUSPARSE_INVALID_VALUE              -3003
#define MAGMA_tally4_ERR_CUSPARSE_ARCH_MISMATCH              -3004
#define MAGMA_tally4_ERR_CUSPARSE_MAPPING_ERROR              -3005
#define MAGMA_tally4_ERR_CUSPARSE_EXECUTION_FAILED           -3006
#define MAGMA_tally4_ERR_CUSPARSE_INTERNAL_ERROR             -3007
#define MAGMA_tally4_ERR_CUSPARSE_MATRIX_TYPE_NOT_SUPPORTED  -3008
#define MAGMA_tally4_ERR_CUSPARSE_ZERO_PIVOT                 -3009


// ----------------------------------------
// parameter constants
// numbering is consistent with CBLAS and PLASMA; see plasma/include/plasma.h
// also with lapack_cwrapper/include/lapack_enum.h
typedef enum {
    Magma_tally4False         = 0,
    Magma_tally4True          = 1
} magma_tally4_bool_t;

typedef enum {
    Magma_tally4RowMajor      = 101,
    Magma_tally4ColMajor      = 102
} magma_tally4_order_t;

// Magma_tally4_ConjTrans is an alias for those rare occasions (zlarfb, zun*, zher*k)
// where we want Magma_tally4_ConjTrans to convert to Magma_tally4Trans in precision generation.
typedef enum {
    Magma_tally4NoTrans       = 111,
    Magma_tally4Trans         = 112,
    Magma_tally4ConjTrans     = 113,
    Magma_tally4_ConjTrans    = Magma_tally4ConjTrans
} magma_tally4_trans_t;

typedef enum {
    Magma_tally4Upper         = 121,
    Magma_tally4Lower         = 122,
    Magma_tally4UpperLower    = 123,
    Magma_tally4Full          = 123,  /* lascl, laset */
    Magma_tally4Hessenberg    = 124   /* lascl */
} magma_tally4_uplo_t;

typedef magma_tally4_uplo_t magma_tally4_type_t;  /* lascl */

typedef enum {
    Magma_tally4NonUnit       = 131,
    Magma_tally4Unit          = 132
} magma_tally4_diag_t;

typedef enum {
    Magma_tally4Left          = 141,
    Magma_tally4Right         = 142,
    Magma_tally4BothSides     = 143   /* trevc */
} magma_tally4_side_t;

typedef enum {
    Magma_tally4OneNorm       = 171,  /* lange, lanhe */
    Magma_tally4RealOneNorm   = 172,
    Magma_tally4TwoNorm       = 173,
    Magma_tally4FrobeniusNorm = 174,
    Magma_tally4InfNorm       = 175,
    Magma_tally4RealInfNorm   = 176,
    Magma_tally4MaxNorm       = 177,
    Magma_tally4RealMaxNorm   = 178
} magma_tally4_norm_t;

typedef enum {
    Magma_tally4DistUniform   = 201,  /* latms */
    Magma_tally4DistSymmetric = 202,
    Magma_tally4DistNormal    = 203
} magma_tally4_dist_t;

typedef enum {
    Magma_tally4HermGeev      = 241,  /* latms */
    Magma_tally4HermPoev      = 242,
    Magma_tally4NonsymPosv    = 243,
    Magma_tally4SymPosv       = 244
} magma_tally4_sym_t;

typedef enum {
    Magma_tally4NoPacking     = 291,  /* latms */
    Magma_tally4PackSubdiag   = 292,
    Magma_tally4PackSupdiag   = 293,
    Magma_tally4PackColumn    = 294,
    Magma_tally4PackRow       = 295,
    Magma_tally4PackLowerBand = 296,
    Magma_tally4PackUpeprBand = 297,
    Magma_tally4PackAll       = 298
} magma_tally4_pack_t;

typedef enum {
    Magma_tally4NoVec         = 301,  /* geev, syev, gesvd */
    Magma_tally4Vec           = 302,  /* geev, syev */
    Magma_tally4IVec          = 303,  /* stedc */
    Magma_tally4AllVec        = 304,  /* gesvd, trevc */
    Magma_tally4SomeVec       = 305,  /* gesvd, trevc */
    Magma_tally4OverwriteVec  = 306,  /* gesvd */
    Magma_tally4BacktransVec  = 307   /* trevc */
} magma_tally4_vec_t;

typedef enum {
    Magma_tally4RangeAll      = 311,  /* syevx, etc. */
    Magma_tally4RangeV        = 312,
    Magma_tally4RangeI        = 313
} magma_tally4_range_t;

typedef enum {
    Magma_tally4Q             = 322,  /* unmbr, ungbr */
    Magma_tally4P             = 323
} magma_tally4_vect_t;

typedef enum {
    Magma_tally4Forward       = 391,  /* larfb */
    Magma_tally4Backward      = 392
} magma_tally4_direct_t;

typedef enum {
    Magma_tally4Columnwise    = 401,  /* larfb */
    Magma_tally4Rowwise       = 402
} magma_tally4_storev_t;

// --------------------
// sparse
typedef enum {
    Magma_tally4_CSR          = 411,
    Magma_tally4_ELLPACKT     = 412,
    Magma_tally4_ELL          = 413,
    Magma_tally4_DENSE        = 414,
    Magma_tally4_BCSR         = 415,
    Magma_tally4_CSC          = 416,
    Magma_tally4_HYB          = 417,
    Magma_tally4_COO          = 418,
    Magma_tally4_ELLRT        = 419,
    Magma_tally4_SPMVFUNCTION = 420,
    Magma_tally4_SELLP        = 421,
    Magma_tally4_ELLD         = 422,

    Magma_tally4_CSRD         = 424,
    Magma_tally4_CSRL         = 427,
    Magma_tally4_CSRU         = 428,
    Magma_tally4_CSRCOO       = 429
} magma_tally4_storage_t;


typedef enum {
    Magma_tally4_CG           = 431,
    Magma_tally4_CGMERGE      = 432,
    Magma_tally4_GMRES        = 433,
    Magma_tally4_BICGSTAB     = 434,
  Magma_tally4_BICGSTABMERGE  = 435,
  Magma_tally4_BICGSTABMERGE2 = 436,
    Magma_tally4_JACOBI       = 437,
    Magma_tally4_GS           = 438,
    Magma_tally4_ITERREF      = 439,
    Magma_tally4_BCSRLU       = 440,
    Magma_tally4_PCG          = 441,
    Magma_tally4_PGMRES       = 442,
    Magma_tally4_PBICGSTAB    = 443,
    Magma_tally4_PASTIX       = 444,
    Magma_tally4_ILU          = 445,
    Magma_tally4_ICC          = 446,
    Magma_tally4_AILU         = 447,
    Magma_tally4_AICC         = 448,
    Magma_tally4_BAITER       = 449,
    Magma_tally4_LOBPCG       = 450,
    Magma_tally4_NONE         = 451,
    Magma_tally4_FUNCTION     = 452
} magma_tally4_solver_type;

typedef enum {
    Magma_tally4_CGS          = 461,
    Magma_tally4_FUSED_CGS    = 462,
    Magma_tally4_MGS          = 463
} magma_tally4_ortho_t;

typedef enum {
    Magma_tally4_CPU          = 471,
    Magma_tally4_DEV          = 472
} magma_tally4_location_t;

typedef enum {
    Magma_tally4_GENERAL      = 481,
    Magma_tally4_SYMMETRIC    = 482
} magma_tally4_symmetry_t;

typedef enum {
    Magma_tally4_ORDERED      = 491,
    Magma_tally4_DIAGFIRST    = 492,
    Magma_tally4_UNITY        = 493,
    Magma_tally4_VALUE        = 494
} magma_tally4_diagorder_t;

typedef enum {
    Magma_tally4_DCOMPLEX     = 501,
    Magma_tally4_FCOMPLEX     = 502,
    Magma_tally4_DOUBLE       = 503,
    Magma_tally4_FLOAT        = 504
} magma_tally4_precision;

typedef enum {
    Magma_tally4_NOSCALE      = 511,
    Magma_tally4_UNITROW      = 512,
    Magma_tally4_UNITDIAG     = 513
} magma_tally4_scale_t;

typedef enum {
    Magma_tally4_FULL         = 521,
    Magma_tally4_LOWER        = 522,
    Magma_tally4_UPPER        = 523
} magma_tally4_fillmode_t;



// When adding constants, remember to do these steps as appropriate:
// 1)  add magma_tally4_xxxx_const()  converter below and in control/constants.cpp
// 2a) add to magma_tally42lapack_const_tally4ants[] in control/constants.cpp
// 2b) update min & max here, which are used to check bounds for magma_tally42lapack_const_tally4ants[]
// 2c) add lapack_xxxx_const() converter below and in control/constants.cpp
#define Magma_tally42lapack_Min  Magma_tally4False     // 0
#define Magma_tally42lapack_Max  Magma_tally4Rowwise   // 402


// ----------------------------------------
// string constants for calling Fortran BLAS and LAPACK
// todo: use translators instead? lapack_const_tally4( Magma_tally4Upper )
#define Magma_tally4RowMajorStr      "Row"
#define Magma_tally4ColMajorStr      "Col"

#define Magma_tally4NoTransStr       "NoTrans"
#define Magma_tally4TransStr         "Trans"
#define Magma_tally4ConjTransStr     "ConjTrans"
#define Magma_tally4_ConjTransStr    "ConjTrans"

#define Magma_tally4UpperStr         "Upper"
#define Magma_tally4LowerStr         "Lower"
#define Magma_tally4UpperLowerStr    "Full"
#define Magma_tally4FullStr          "Full"

#define Magma_tally4NonUnitStr       "NonUnit"
#define Magma_tally4UnitStr          "Unit"

#define Magma_tally4LeftStr          "Left"
#define Magma_tally4RightStr         "Right"
#define Magma_tally4BothSidesStr     "Both"

#define Magma_tally4OneNormStr       "1"
#define Magma_tally4TwoNormStr       "2"
#define Magma_tally4FrobeniusNormStr "Fro"
#define Magma_tally4InfNormStr       "Inf"
#define Magma_tally4MaxNormStr       "Max"

#define Magma_tally4ForwardStr       "Forward"
#define Magma_tally4BackwardStr      "Backward"

#define Magma_tally4ColumnwiseStr    "Columnwise"
#define Magma_tally4RowwiseStr       "Rowwise"

#define Magma_tally4NoVecStr         "NoVec"
#define Magma_tally4VecStr           "Vec"
#define Magma_tally4IVecStr          "IVec"
#define Magma_tally4AllVecStr        "All"
#define Magma_tally4SomeVecStr       "Some"
#define Magma_tally4OverwriteVecStr  "Overwrite"


#ifdef __cplusplus
extern "C" {
#endif

// --------------------
// Convert LAPACK character constants to MAGMA_tally4 constants.
// This is a one-to-many mapping, requiring multiple translators
// (e.g., "N" can be NoTrans or NonUnit or NoVec).
magma_tally4_bool_t   magma_tally4_bool_const  ( char lapack_char );
magma_tally4_order_t  magma_tally4_order_const ( char lapack_char );
magma_tally4_trans_t  magma_tally4_trans_const ( char lapack_char );
magma_tally4_uplo_t   magma_tally4_uplo_const  ( char lapack_char );
magma_tally4_diag_t   magma_tally4_diag_const  ( char lapack_char );
magma_tally4_side_t   magma_tally4_side_const  ( char lapack_char );
magma_tally4_norm_t   magma_tally4_norm_const  ( char lapack_char );
magma_tally4_dist_t   magma_tally4_dist_const  ( char lapack_char );
magma_tally4_sym_t    magma_tally4_sym_const   ( char lapack_char );
magma_tally4_pack_t   magma_tally4_pack_const  ( char lapack_char );
magma_tally4_vec_t    magma_tally4_vec_const   ( char lapack_char );
magma_tally4_range_t  magma_tally4_range_const ( char lapack_char );
magma_tally4_vect_t   magma_tally4_vect_const  ( char lapack_char );
magma_tally4_direct_t magma_tally4_direct_const( char lapack_char );
magma_tally4_storev_t magma_tally4_storev_const( char lapack_char );


// --------------------
// Convert MAGMA_tally4 constants to LAPACK(E) constants.
// The generic lapack_const_tally4 works for all cases, but the specific routines
// (e.g., lapack_trans_const_tally4) do better error checking.
const char* lapack_const_tally4       ( int            magma_tally4_const );
const char* lapack_bool_const_tally4  ( magma_tally4_bool_t   magma_tally4_const );
const char* lapack_order_const_tally4 ( magma_tally4_order_t  magma_tally4_const );
const char* lapack_trans_const_tally4 ( magma_tally4_trans_t  magma_tally4_const );
const char* lapack_uplo_const_tally4  ( magma_tally4_uplo_t   magma_tally4_const );
const char* lapack_diag_const_tally4  ( magma_tally4_diag_t   magma_tally4_const );
const char* lapack_side_const_tally4  ( magma_tally4_side_t   magma_tally4_const );
const char* lapack_norm_const_tally4  ( magma_tally4_norm_t   magma_tally4_const );
const char* lapack_dist_const_tally4  ( magma_tally4_dist_t   magma_tally4_const );
const char* lapack_sym_const_tally4   ( magma_tally4_sym_t    magma_tally4_const );
const char* lapack_pack_const_tally4  ( magma_tally4_pack_t   magma_tally4_const );
const char* lapack_vec_const_tally4   ( magma_tally4_vec_t    magma_tally4_const );
const char* lapack_range_const_tally4 ( magma_tally4_range_t  magma_tally4_const );
const char* lapack_vect_const_tally4  ( magma_tally4_vect_t   magma_tally4_const );
const char* lapack_direct_const_tally4( magma_tally4_direct_t magma_tally4_const );
const char* lapack_storev_const_tally4( magma_tally4_storev_t magma_tally4_const );

static inline char lapacke_const_tally4       ( int magma_tally4_const            ) { return *lapack_const_tally4       ( magma_tally4_const ); }
static inline char lapacke_bool_const_tally4  ( magma_tally4_bool_t   magma_tally4_const ) { return *lapack_bool_const_tally4  ( magma_tally4_const ); }
static inline char lapacke_order_const_tally4 ( magma_tally4_order_t  magma_tally4_const ) { return *lapack_order_const_tally4 ( magma_tally4_const ); }
static inline char lapacke_trans_const_tally4 ( magma_tally4_trans_t  magma_tally4_const ) { return *lapack_trans_const_tally4 ( magma_tally4_const ); }
static inline char lapacke_uplo_const_tally4  ( magma_tally4_uplo_t   magma_tally4_const ) { return *lapack_uplo_const_tally4  ( magma_tally4_const ); }
static inline char lapacke_diag_const_tally4  ( magma_tally4_diag_t   magma_tally4_const ) { return *lapack_diag_const_tally4  ( magma_tally4_const ); }
static inline char lapacke_side_const_tally4  ( magma_tally4_side_t   magma_tally4_const ) { return *lapack_side_const_tally4  ( magma_tally4_const ); }
static inline char lapacke_norm_const_tally4  ( magma_tally4_norm_t   magma_tally4_const ) { return *lapack_norm_const_tally4  ( magma_tally4_const ); }
static inline char lapacke_dist_const_tally4  ( magma_tally4_dist_t   magma_tally4_const ) { return *lapack_dist_const_tally4  ( magma_tally4_const ); }
static inline char lapacke_sym_const_tally4   ( magma_tally4_sym_t    magma_tally4_const ) { return *lapack_sym_const_tally4   ( magma_tally4_const ); }
static inline char lapacke_pack_const_tally4  ( magma_tally4_pack_t   magma_tally4_const ) { return *lapack_pack_const_tally4  ( magma_tally4_const ); }
static inline char lapacke_vec_const_tally4   ( magma_tally4_vec_t    magma_tally4_const ) { return *lapack_vec_const_tally4   ( magma_tally4_const ); }
static inline char lapacke_range_const_tally4 ( magma_tally4_range_t  magma_tally4_const ) { return *lapack_range_const_tally4 ( magma_tally4_const ); }
static inline char lapacke_vect_const_tally4  ( magma_tally4_vect_t   magma_tally4_const ) { return *lapack_vect_const_tally4  ( magma_tally4_const ); }
static inline char lapacke_direct_const_tally4( magma_tally4_direct_t magma_tally4_const ) { return *lapack_direct_const_tally4( magma_tally4_const ); }
static inline char lapacke_storev_const_tally4( magma_tally4_storev_t magma_tally4_const ) { return *lapack_storev_const_tally4( magma_tally4_const ); }


// --------------------
// Convert MAGMA_tally4 constants to clAmdBlas constants.
#if defined(HAVE_clAmdBlas)
clAmdBlasOrder       amdblas_order_const( magma_tally4_order_t order );
clAmdBlasTranspose   amdblas_trans_const( magma_tally4_trans_t trans );
clAmdBlasUplo        amdblas_uplo_const ( magma_tally4_uplo_t  uplo  );
clAmdBlasDiag        amdblas_diag_const ( magma_tally4_diag_t  diag  );
clAmdBlasSide        amdblas_side_const ( magma_tally4_side_t  side  );
#endif


// --------------------
// Convert MAGMA_tally4 constants to CUBLAS constants.
#if defined(CUBLAS_V2_H_)
cublasOperation_t    cublas_trans_const_tally4 ( magma_tally4_trans_t trans );
cublasFillMode_t     cublas_uplo_const_tally4  ( magma_tally4_uplo_t  uplo  );
cublasDiagType_t     cublas_diag_const_tally4  ( magma_tally4_diag_t  diag  );
cublasSideMode_t     cublas_side_const_tally4  ( magma_tally4_side_t  side  );
#endif


// --------------------
// Convert MAGMA_tally4 constants to CBLAS constants.
#if defined(HAVE_CBLAS)
#include <cblas.h>
enum CBLAS_ORDER     cblas_order_const  ( magma_tally4_order_t order );
enum CBLAS_TRANSPOSE cblas_trans_const  ( magma_tally4_trans_t trans );
enum CBLAS_UPLO      cblas_uplo_const   ( magma_tally4_uplo_t  uplo  );
enum CBLAS_DIAG      cblas_diag_const   ( magma_tally4_diag_t  diag  );
enum CBLAS_SIDE      cblas_side_const   ( magma_tally4_side_t  side  );
#endif


#ifdef __cplusplus
}
#endif

#endif        //  #ifndef MAGMA_tally4_TYPES_H
