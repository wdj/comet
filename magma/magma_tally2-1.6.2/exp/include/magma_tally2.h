/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
*/

#include <quark.h>

#ifndef _MAGMA_tally2_
#define _MAGMA_tally2_

/* ------------------------------------------------------------
 * MAGMA_tally2 Blas Functions 
 * --------------------------------------------------------- */ 
#include "magma_tally2blas.h"

#include "auxiliary.h"

/* ------------------------------------------------------------
 * MAGMA_tally2 Context
 * --------------------------------------------------------- */

typedef struct context
{
  /* Number of CPU core in this context */
  magma_tally2_int_t num_cores;

  /* Number of GPUs in this context */
  magma_tally2_int_t num_gpus;

  /* GPU contexts */
  CUcontext *gpu_context;

  /* QUARK scheduler */
  Quark *quark;

  /* Block size, internally used for some algorithms */
  magma_tally2_int_t nb;

  /* Pointer to other global algorithm-dependent parameters */
  void *params;

} magma_tally2_context;

/* ------------------------------------------------------------
 * MAGMA_tally2 functions
 * --------------------------------------------------------- */
#include "magma_tally2_z.h"
#include "magma_tally2_c.h"
#include "magma_tally2_d.h"
#include "magma_tally2_s.h"
#include "magma_tally2_zc.h"
#include "magma_tally2_ds.h"

#define Magma_tally2NoTrans       'N'
#define Magma_tally2Trans         'T'
#define Magma_tally2ConjTrans     'C'

#define Magma_tally2Upper         'U'
#define Magma_tally2Lower         'L'
#define Magma_tally2UpperLower    'A'

#define Magma_tally2NonUnit       'N'
#define Magma_tally2Unit          'U'

#define Magma_tally2Left          'L'
#define Magma_tally2Right         'R'

#define Magma_tally2Forward       'F'
#define Magma_tally2Backward      'B'
                           
#define Magma_tally2Columnwise    'C'
#define Magma_tally2Rowwise       'R'

#define Magma_tally2NoVectors     'N'
#define Magma_tally2Vectors       'V'

#define Magma_tally2NoTransStr    "NonTrans"
#define Magma_tally2TransStr      "Trans"
#define Magma_tally2ConjTransStr  "Conj"

#define Magma_tally2UpperStr      "Upper"
#define Magma_tally2LowerStr      "Lower"
#define Magma_tally2UpperLowerStr "All"

#define Magma_tally2NonUnitStr    "NonUnit"
#define Magma_tally2UnitStr       "Unit"

#define Magma_tally2LeftStr       "Left"
#define Magma_tally2RightStr      "Right"

#define Magma_tally2ForwardStr    "Forward"
#define Magma_tally2BackwardStr   "Backward"

#define Magma_tally2ColumnwiseStr "Columnwise"
#define Magma_tally2RowwiseStr    "Rowwise"

#define Magma_tally2NoVectorsStr  "NoVectors"
#define Magma_tally2VectorsStr    "Vectors"

/* ------------------------------------------------------------
 *   Return codes
 * --------------------------------------------------------- */
#define MAGMA_tally2_SUCCESS             0
#define MAGMA_tally2_ERR_ILLEGAL_VALUE  -4
#define MAGMA_tally2_ERR_ALLOCATION     -5
#define MAGMA_tally2_ERR_HOSTALLOC      -6
#define MAGMA_tally2_ERR_CUBLASALLOC    -7

/* ------------------------------------------------------------
 *   Macros to deal with cuda complex
 * --------------------------------------------------------- */
#define MAGMA_tally2_Z_SET2REAL(v, t)    (v).x = (t); (v).y = 0.0
#define MAGMA_tally2_Z_OP_NEG_ASGN(t, z) (t).x = -(z).x; (t).y = -(z).y
#define MAGMA_tally2_Z_EQUAL(u,v)        (((u).x == (v).x) && ((u).y == (v).y))
#define MAGMA_tally2_Z_GET_X(u)          ((u).x)
#define MAGMA_tally2_Z_ASSIGN(v, t)      (v).x = (t).x; (v).y = (t).y
#define MAGMA_tally2_Z_CNJG(v, t)        (v).x = (t).x; (v).y = -(t).y
#define MAGMA_tally2_Z_DSCALE(v, t, s)   (v).x = (t).x/(s); (v).y = (t).y/(s)      
#define MAGMA_tally2_Z_OP_NEG(a, b, c)   (a).x = (b).x-(c).x; (a).y = (b).y-(c).y
#define MAGMA_tally2_Z_MAKE(r, i)        make_cuDoubleComplex((r), (i))
#define MAGMA_tally2_Z_REAL(a)           cuCreal(a)
#define MAGMA_tally2_Z_IMAG(a)           cuCimag(a)
#define MAGMA_tally2_Z_ADD(a, b)         cuCadd((a), (b))
#define MAGMA_tally2_Z_SUB(a, b)         cuCsub((a), (b))
#define MAGMA_tally2_Z_MUL(a, b)         cuCmul((a), (b))
#define MAGMA_tally2_Z_DIV(a, b)         cuCdiv((a), (b))
#define MAGMA_tally2_Z_ABS(a)            cuCabs((a))
#define MAGMA_tally2_Z_ZERO              make_cuDoubleComplex(0.0, 0.0)
#define MAGMA_tally2_Z_ONE               make_cuDoubleComplex(1.0, 0.0)
#define MAGMA_tally2_Z_HALF              make_cuDoubleComplex(0.5, 0.0)
#define MAGMA_tally2_Z_NEG_ONE           make_cuDoubleComplex(-1.0, 0.0)
#define MAGMA_tally2_Z_NEG_HALF          make_cuDoubleComplex(-0.5, 0.0)

#define MAGMA_tally2_C_SET2REAL(v, t)    (v).x = (t); (v).y = 0.0
#define MAGMA_tally2_C_OP_NEG_ASGN(t, z) (t).x = -(z).x; (t).y = -(z).y
#define MAGMA_tally2_C_EQUAL(u,v)        (((u).x == (v).x) && ((u).y == (v).y))
#define MAGMA_tally2_C_GET_X(u)          ((u).x)
#define MAGMA_tally2_C_ASSIGN(v, t)      (v).x = (t).x; (v).y = (t).y
#define MAGMA_tally2_C_CNJG(v, t)        (v).x= (t).x; (v).y = -(t).y
#define MAGMA_tally2_C_DSCALE(v, t, s)   (v).x = (t).x/(s); (v).y = (t).y/(s)
#define MAGMA_tally2_C_OP_NEG(a, b, c)   (a).x = (b).x-(c).x; (a).y = (b).y-(c).y
#define MAGMA_tally2_C_MAKE(r, i)        make_cuFloatComplex((r), (i))
#define MAGMA_tally2_C_REAL(a)           cuCrealf(a)
#define MAGMA_tally2_C_IMAG(a)           cuCimagf(a)
#define MAGMA_tally2_C_ADD(a, b)         cuCaddf((a), (b))
#define MAGMA_tally2_C_SUB(a, b)         cuCsubf((a), (b))
#define MAGMA_tally2_C_MUL(a, b)         cuCmulf((a), (b))
#define MAGMA_tally2_C_DIV(a, b)         cuCdivf((a), (b))
#define MAGMA_tally2_C_ABS(a)            cuCabsf((a))
#define MAGMA_tally2_C_ZERO              make_cuFloatComplex(0.0, 0.0)
#define MAGMA_tally2_C_ONE               make_cuFloatComplex(1.0, 0.0)
#define MAGMA_tally2_C_HALF              make_cuFloatComplex(0.5, 0.0)
#define MAGMA_tally2_C_NEG_ONE           make_cuFloatComplex(-1.0, 0.0)
#define MAGMA_tally2_C_NEG_HALF          make_cuFloatComplex(-0.5, 0.0)

#define MAGMA_tally2_D_SET2REAL(v, t)    (v) = (t);
#define MAGMA_tally2_D_OP_NEG_ASGN(t, z) (t) = -(z)
#define MAGMA_tally2_D_EQUAL(u,v)        ((u) == (v))
#define MAGMA_tally2_D_GET_X(u)          (u)
#define MAGMA_tally2_D_ASSIGN(v, t)      (v) = (t)
#define MAGMA_tally2_D_CNJG(v, t)        (v) = (t)
#define MAGMA_tally2_D_DSCALE(v, t, s)   (v) = (t)/(s)
#define MAGMA_tally2_D_OP_NEG(a, b, c)   (a) = (b) - (c)
#define MAGMA_tally2_D_MAKE(r, i)        (r)
#define MAGMA_tally2_D_REAL(a)           (a)
#define MAGMA_tally2_D_IMAG(a)           (a)
#define MAGMA_tally2_D_ADD(a, b)         ( (a) + (b) )
#define MAGMA_tally2_D_SUB(a, b)         ( (a) - (b) )
#define MAGMA_tally2_D_MUL(a, b)         ( (a) * (b) )
#define MAGMA_tally2_D_DIV(a, b)         ( (a) / (b) )
#define MAGMA_tally2_D_ABS(a)            ((a)>0?(a):-(a))
#define MAGMA_tally2_D_ZERO              (0.0)
#define MAGMA_tally2_D_ONE               (1.0)
#define MAGMA_tally2_D_HALF              (0.5)
#define MAGMA_tally2_D_NEG_ONE           (-1.0)
#define MAGMA_tally2_D_NEG_HALF          (-0.5)

#define MAGMA_tally2_S_SET2REAL(v, t)    (v) = (t);
#define MAGMA_tally2_S_OP_NEG_ASGN(t, z) (t) = -(z)
#define MAGMA_tally2_S_EQUAL(u,v)        ((u) == (v))
#define MAGMA_tally2_S_GET_X(u)          (u)
#define MAGMA_tally2_S_ASSIGN(v, t)      (v) = (t)
#define MAGMA_tally2_S_CNJG(v, t)        (v) = (t)
#define MAGMA_tally2_S_DSCALE(v, t, s)   (v) = (t)/(s)
#define MAGMA_tally2_S_OP_NEG(a, b, c)   (a) = (b) - (c)
#define MAGMA_tally2_S_MAKE(r, i)        (r)
#define MAGMA_tally2_S_REAL(a)           (a)
#define MAGMA_tally2_S_IMAG(a)           (a)
#define MAGMA_tally2_S_ADD(a, b)         ( (a) + (b) )
#define MAGMA_tally2_S_SUB(a, b)         ( (a) - (b) )
#define MAGMA_tally2_S_MUL(a, b)         ( (a) * (b) )
#define MAGMA_tally2_S_DIV(a, b)         ( (a) / (b) )
#define MAGMA_tally2_S_ABS(a)            ((a)>0?(a):-(a))
#define MAGMA_tally2_S_ZERO              (0.0)
#define MAGMA_tally2_S_ONE               (1.0)
#define MAGMA_tally2_S_HALF              (0.5)
#define MAGMA_tally2_S_NEG_ONE           (-1.0)
#define MAGMA_tally2_S_NEG_HALF          (-0.5)

#ifndef CBLAS_SADDR
#define CBLAS_SADDR(a)  &(a)
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------
 *   -- MAGMA_tally2 function definitions
 * --------------------------------------------------------- */
void magma_tally2_xerbla( const char *name, magma_tally2_int_t info );
magma_tally2_context *magma_tally2_init(void *, void* (*func)(void *a), magma_tally2_int_t nthread, magma_tally2_int_t ncpu, 
                          magma_tally2_int_t ngpu, magma_tally2_int_t argc, char **argv);
void magma_tally2_finalize(magma_tally2_context *cntxt);
void auto_tune(char algorithm, char precision, magma_tally2_int_t ncores, magma_tally2_int_t ncorespsocket,
               magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t *nb, magma_tally2_int_t *ob, magma_tally2_int_t *ib,
               magma_tally2_int_t *nthreads, magma_tally2_int_t *nquarkthreads);



#ifdef __cplusplus
}
#endif

#endif

