/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
*/

#include <quark.h>

#ifndef _MAGMA_minproduct_
#define _MAGMA_minproduct_

/* ------------------------------------------------------------
 * MAGMA_minproduct Blas Functions 
 * --------------------------------------------------------- */ 
#include "magma_minproductblas.h"

#include "auxiliary.h"

/* ------------------------------------------------------------
 * MAGMA_minproduct Context
 * --------------------------------------------------------- */

typedef struct context
{
  /* Number of CPU core in this context */
  magma_minproduct_int_t num_cores;

  /* Number of GPUs in this context */
  magma_minproduct_int_t num_gpus;

  /* GPU contexts */
  CUcontext *gpu_context;

  /* QUARK scheduler */
  Quark *quark;

  /* Block size, internally used for some algorithms */
  magma_minproduct_int_t nb;

  /* Pointer to other global algorithm-dependent parameters */
  void *params;

} magma_minproduct_context;

/* ------------------------------------------------------------
 * MAGMA_minproduct functions
 * --------------------------------------------------------- */
#include "magma_minproduct_z.h"
#include "magma_minproduct_c.h"
#include "magma_minproduct_d.h"
#include "magma_minproduct_s.h"
#include "magma_minproduct_zc.h"
#include "magma_minproduct_ds.h"

#define Magma_minproductNoTrans       'N'
#define Magma_minproductTrans         'T'
#define Magma_minproductConjTrans     'C'

#define Magma_minproductUpper         'U'
#define Magma_minproductLower         'L'
#define Magma_minproductUpperLower    'A'

#define Magma_minproductNonUnit       'N'
#define Magma_minproductUnit          'U'

#define Magma_minproductLeft          'L'
#define Magma_minproductRight         'R'

#define Magma_minproductForward       'F'
#define Magma_minproductBackward      'B'
                           
#define Magma_minproductColumnwise    'C'
#define Magma_minproductRowwise       'R'

#define Magma_minproductNoVectors     'N'
#define Magma_minproductVectors       'V'

#define Magma_minproductNoTransStr    "NonTrans"
#define Magma_minproductTransStr      "Trans"
#define Magma_minproductConjTransStr  "Conj"

#define Magma_minproductUpperStr      "Upper"
#define Magma_minproductLowerStr      "Lower"
#define Magma_minproductUpperLowerStr "All"

#define Magma_minproductNonUnitStr    "NonUnit"
#define Magma_minproductUnitStr       "Unit"

#define Magma_minproductLeftStr       "Left"
#define Magma_minproductRightStr      "Right"

#define Magma_minproductForwardStr    "Forward"
#define Magma_minproductBackwardStr   "Backward"

#define Magma_minproductColumnwiseStr "Columnwise"
#define Magma_minproductRowwiseStr    "Rowwise"

#define Magma_minproductNoVectorsStr  "NoVectors"
#define Magma_minproductVectorsStr    "Vectors"

/* ------------------------------------------------------------
 *   Return codes
 * --------------------------------------------------------- */
#define MAGMA_minproduct_SUCCESS             0
#define MAGMA_minproduct_ERR_ILLEGAL_VALUE  -4
#define MAGMA_minproduct_ERR_ALLOCATION     -5
#define MAGMA_minproduct_ERR_HOSTALLOC      -6
#define MAGMA_minproduct_ERR_CUBLASALLOC    -7

/* ------------------------------------------------------------
 *   Macros to deal with cuda complex
 * --------------------------------------------------------- */
#define MAGMA_minproduct_Z_SET2REAL(v, t)    (v).x = (t); (v).y = 0.0
#define MAGMA_minproduct_Z_OP_NEG_ASGN(t, z) (t).x = -(z).x; (t).y = -(z).y
#define MAGMA_minproduct_Z_EQUAL(u,v)        (((u).x == (v).x) && ((u).y == (v).y))
#define MAGMA_minproduct_Z_GET_X(u)          ((u).x)
#define MAGMA_minproduct_Z_ASSIGN(v, t)      (v).x = (t).x; (v).y = (t).y
#define MAGMA_minproduct_Z_CNJG(v, t)        (v).x = (t).x; (v).y = -(t).y
#define MAGMA_minproduct_Z_DSCALE(v, t, s)   (v).x = (t).x/(s); (v).y = (t).y/(s)      
#define MAGMA_minproduct_Z_OP_NEG(a, b, c)   (a).x = (b).x-(c).x; (a).y = (b).y-(c).y
#define MAGMA_minproduct_Z_MAKE(r, i)        make_cuDoubleComplex((r), (i))
#define MAGMA_minproduct_Z_REAL(a)           cuCreal(a)
#define MAGMA_minproduct_Z_IMAG(a)           cuCimag(a)
#define MAGMA_minproduct_Z_ADD(a, b)         cuCadd((a), (b))
#define MAGMA_minproduct_Z_SUB(a, b)         cuCsub((a), (b))
#define MAGMA_minproduct_Z_MUL(a, b)         cuCmul((a), (b))
#define MAGMA_minproduct_Z_DIV(a, b)         cuCdiv((a), (b))
#define MAGMA_minproduct_Z_ABS(a)            cuCabs((a))
#define MAGMA_minproduct_Z_ZERO              make_cuDoubleComplex(0.0, 0.0)
#define MAGMA_minproduct_Z_ONE               make_cuDoubleComplex(1.0, 0.0)
#define MAGMA_minproduct_Z_HALF              make_cuDoubleComplex(0.5, 0.0)
#define MAGMA_minproduct_Z_NEG_ONE           make_cuDoubleComplex(-1.0, 0.0)
#define MAGMA_minproduct_Z_NEG_HALF          make_cuDoubleComplex(-0.5, 0.0)

#define MAGMA_minproduct_C_SET2REAL(v, t)    (v).x = (t); (v).y = 0.0
#define MAGMA_minproduct_C_OP_NEG_ASGN(t, z) (t).x = -(z).x; (t).y = -(z).y
#define MAGMA_minproduct_C_EQUAL(u,v)        (((u).x == (v).x) && ((u).y == (v).y))
#define MAGMA_minproduct_C_GET_X(u)          ((u).x)
#define MAGMA_minproduct_C_ASSIGN(v, t)      (v).x = (t).x; (v).y = (t).y
#define MAGMA_minproduct_C_CNJG(v, t)        (v).x= (t).x; (v).y = -(t).y
#define MAGMA_minproduct_C_DSCALE(v, t, s)   (v).x = (t).x/(s); (v).y = (t).y/(s)
#define MAGMA_minproduct_C_OP_NEG(a, b, c)   (a).x = (b).x-(c).x; (a).y = (b).y-(c).y
#define MAGMA_minproduct_C_MAKE(r, i)        make_cuFloatComplex((r), (i))
#define MAGMA_minproduct_C_REAL(a)           cuCrealf(a)
#define MAGMA_minproduct_C_IMAG(a)           cuCimagf(a)
#define MAGMA_minproduct_C_ADD(a, b)         cuCaddf((a), (b))
#define MAGMA_minproduct_C_SUB(a, b)         cuCsubf((a), (b))
#define MAGMA_minproduct_C_MUL(a, b)         cuCmulf((a), (b))
#define MAGMA_minproduct_C_DIV(a, b)         cuCdivf((a), (b))
#define MAGMA_minproduct_C_ABS(a)            cuCabsf((a))
#define MAGMA_minproduct_C_ZERO              make_cuFloatComplex(0.0, 0.0)
#define MAGMA_minproduct_C_ONE               make_cuFloatComplex(1.0, 0.0)
#define MAGMA_minproduct_C_HALF              make_cuFloatComplex(0.5, 0.0)
#define MAGMA_minproduct_C_NEG_ONE           make_cuFloatComplex(-1.0, 0.0)
#define MAGMA_minproduct_C_NEG_HALF          make_cuFloatComplex(-0.5, 0.0)

#define MAGMA_minproduct_D_SET2REAL(v, t)    (v) = (t);
#define MAGMA_minproduct_D_OP_NEG_ASGN(t, z) (t) = -(z)
#define MAGMA_minproduct_D_EQUAL(u,v)        ((u) == (v))
#define MAGMA_minproduct_D_GET_X(u)          (u)
#define MAGMA_minproduct_D_ASSIGN(v, t)      (v) = (t)
#define MAGMA_minproduct_D_CNJG(v, t)        (v) = (t)
#define MAGMA_minproduct_D_DSCALE(v, t, s)   (v) = (t)/(s)
#define MAGMA_minproduct_D_OP_NEG(a, b, c)   (a) = (b) - (c)
#define MAGMA_minproduct_D_MAKE(r, i)        (r)
#define MAGMA_minproduct_D_REAL(a)           (a)
#define MAGMA_minproduct_D_IMAG(a)           (a)
#define MAGMA_minproduct_D_ADD(a, b)         ( (a) + (b) )
#define MAGMA_minproduct_D_SUB(a, b)         ( (a) - (b) )
#define MAGMA_minproduct_D_MUL(a, b)         ( (a) * (b) )
#define MAGMA_minproduct_D_DIV(a, b)         ( (a) / (b) )
#define MAGMA_minproduct_D_ABS(a)            ((a)>0?(a):-(a))
#define MAGMA_minproduct_D_ZERO              (0.0)
#define MAGMA_minproduct_D_ONE               (1.0)
#define MAGMA_minproduct_D_HALF              (0.5)
#define MAGMA_minproduct_D_NEG_ONE           (-1.0)
#define MAGMA_minproduct_D_NEG_HALF          (-0.5)

#define MAGMA_minproduct_S_SET2REAL(v, t)    (v) = (t);
#define MAGMA_minproduct_S_OP_NEG_ASGN(t, z) (t) = -(z)
#define MAGMA_minproduct_S_EQUAL(u,v)        ((u) == (v))
#define MAGMA_minproduct_S_GET_X(u)          (u)
#define MAGMA_minproduct_S_ASSIGN(v, t)      (v) = (t)
#define MAGMA_minproduct_S_CNJG(v, t)        (v) = (t)
#define MAGMA_minproduct_S_DSCALE(v, t, s)   (v) = (t)/(s)
#define MAGMA_minproduct_S_OP_NEG(a, b, c)   (a) = (b) - (c)
#define MAGMA_minproduct_S_MAKE(r, i)        (r)
#define MAGMA_minproduct_S_REAL(a)           (a)
#define MAGMA_minproduct_S_IMAG(a)           (a)
#define MAGMA_minproduct_S_ADD(a, b)         ( (a) + (b) )
#define MAGMA_minproduct_S_SUB(a, b)         ( (a) - (b) )
#define MAGMA_minproduct_S_MUL(a, b)         ( (a) * (b) )
#define MAGMA_minproduct_S_DIV(a, b)         ( (a) / (b) )
#define MAGMA_minproduct_S_ABS(a)            ((a)>0?(a):-(a))
#define MAGMA_minproduct_S_ZERO              (0.0)
#define MAGMA_minproduct_S_ONE               (1.0)
#define MAGMA_minproduct_S_HALF              (0.5)
#define MAGMA_minproduct_S_NEG_ONE           (-1.0)
#define MAGMA_minproduct_S_NEG_HALF          (-0.5)

#ifndef CBLAS_SADDR
#define CBLAS_SADDR(a)  &(a)
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------
 *   -- MAGMA_minproduct function definitions
 * --------------------------------------------------------- */
void magma_minproduct_xerbla( const char *name, magma_minproduct_int_t info );
magma_minproduct_context *magma_minproduct_init(void *, void* (*func)(void *a), magma_minproduct_int_t nthread, magma_minproduct_int_t ncpu, 
                          magma_minproduct_int_t ngpu, magma_minproduct_int_t argc, char **argv);
void magma_minproduct_finalize(magma_minproduct_context *cntxt);
void auto_tune(char algorithm, char precision, magma_minproduct_int_t ncores, magma_minproduct_int_t ncorespsocket,
               magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t *nb, magma_minproduct_int_t *ob, magma_minproduct_int_t *ib,
               magma_minproduct_int_t *nthreads, magma_minproduct_int_t *nquarkthreads);



#ifdef __cplusplus
}
#endif

#endif

