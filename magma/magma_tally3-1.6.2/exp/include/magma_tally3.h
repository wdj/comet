/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
*/

#include <quark.h>

#ifndef _MAGMA_tally3_
#define _MAGMA_tally3_

/* ------------------------------------------------------------
 * MAGMA_tally3 Blas Functions 
 * --------------------------------------------------------- */ 
#include "magma_tally3blas.h"

#include "auxiliary.h"

/* ------------------------------------------------------------
 * MAGMA_tally3 Context
 * --------------------------------------------------------- */

typedef struct context
{
  /* Number of CPU core in this context */
  magma_tally3_int_t num_cores;

  /* Number of GPUs in this context */
  magma_tally3_int_t num_gpus;

  /* GPU contexts */
  CUcontext *gpu_context;

  /* QUARK scheduler */
  Quark *quark;

  /* Block size, internally used for some algorithms */
  magma_tally3_int_t nb;

  /* Pointer to other global algorithm-dependent parameters */
  void *params;

} magma_tally3_context;

/* ------------------------------------------------------------
 * MAGMA_tally3 functions
 * --------------------------------------------------------- */
#include "magma_tally3_z.h"
#include "magma_tally3_c.h"
#include "magma_tally3_d.h"
#include "magma_tally3_s.h"
#include "magma_tally3_zc.h"
#include "magma_tally3_ds.h"

#define Magma_tally3NoTrans       'N'
#define Magma_tally3Trans         'T'
#define Magma_tally3ConjTrans     'C'

#define Magma_tally3Upper         'U'
#define Magma_tally3Lower         'L'
#define Magma_tally3UpperLower    'A'

#define Magma_tally3NonUnit       'N'
#define Magma_tally3Unit          'U'

#define Magma_tally3Left          'L'
#define Magma_tally3Right         'R'

#define Magma_tally3Forward       'F'
#define Magma_tally3Backward      'B'
                           
#define Magma_tally3Columnwise    'C'
#define Magma_tally3Rowwise       'R'

#define Magma_tally3NoVectors     'N'
#define Magma_tally3Vectors       'V'

#define Magma_tally3NoTransStr    "NonTrans"
#define Magma_tally3TransStr      "Trans"
#define Magma_tally3ConjTransStr  "Conj"

#define Magma_tally3UpperStr      "Upper"
#define Magma_tally3LowerStr      "Lower"
#define Magma_tally3UpperLowerStr "All"

#define Magma_tally3NonUnitStr    "NonUnit"
#define Magma_tally3UnitStr       "Unit"

#define Magma_tally3LeftStr       "Left"
#define Magma_tally3RightStr      "Right"

#define Magma_tally3ForwardStr    "Forward"
#define Magma_tally3BackwardStr   "Backward"

#define Magma_tally3ColumnwiseStr "Columnwise"
#define Magma_tally3RowwiseStr    "Rowwise"

#define Magma_tally3NoVectorsStr  "NoVectors"
#define Magma_tally3VectorsStr    "Vectors"

/* ------------------------------------------------------------
 *   Return codes
 * --------------------------------------------------------- */
#define MAGMA_tally3_SUCCESS             0
#define MAGMA_tally3_ERR_ILLEGAL_VALUE  -4
#define MAGMA_tally3_ERR_ALLOCATION     -5
#define MAGMA_tally3_ERR_HOSTALLOC      -6
#define MAGMA_tally3_ERR_CUBLASALLOC    -7

/* ------------------------------------------------------------
 *   Macros to deal with cuda complex
 * --------------------------------------------------------- */
#define MAGMA_tally3_Z_SET2REAL(v, t)    (v).x = (t); (v).y = 0.0
#define MAGMA_tally3_Z_OP_NEG_ASGN(t, z) (t).x = -(z).x; (t).y = -(z).y
#define MAGMA_tally3_Z_EQUAL(u,v)        (((u).x == (v).x) && ((u).y == (v).y))
#define MAGMA_tally3_Z_GET_X(u)          ((u).x)
#define MAGMA_tally3_Z_ASSIGN(v, t)      (v).x = (t).x; (v).y = (t).y
#define MAGMA_tally3_Z_CNJG(v, t)        (v).x = (t).x; (v).y = -(t).y
#define MAGMA_tally3_Z_DSCALE(v, t, s)   (v).x = (t).x/(s); (v).y = (t).y/(s)      
#define MAGMA_tally3_Z_OP_NEG(a, b, c)   (a).x = (b).x-(c).x; (a).y = (b).y-(c).y
#define MAGMA_tally3_Z_MAKE(r, i)        make_cuDoubleComplex((r), (i))
#define MAGMA_tally3_Z_REAL(a)           cuCreal(a)
#define MAGMA_tally3_Z_IMAG(a)           cuCimag(a)
#define MAGMA_tally3_Z_ADD(a, b)         cuCadd((a), (b))
#define MAGMA_tally3_Z_SUB(a, b)         cuCsub((a), (b))
#define MAGMA_tally3_Z_MUL(a, b)         cuCmul((a), (b))
#define MAGMA_tally3_Z_DIV(a, b)         cuCdiv((a), (b))
#define MAGMA_tally3_Z_ABS(a)            cuCabs((a))
#define MAGMA_tally3_Z_ZERO              make_cuDoubleComplex(0.0, 0.0)
#define MAGMA_tally3_Z_ONE               make_cuDoubleComplex(1.0, 0.0)
#define MAGMA_tally3_Z_HALF              make_cuDoubleComplex(0.5, 0.0)
#define MAGMA_tally3_Z_NEG_ONE           make_cuDoubleComplex(-1.0, 0.0)
#define MAGMA_tally3_Z_NEG_HALF          make_cuDoubleComplex(-0.5, 0.0)

#define MAGMA_tally3_C_SET2REAL(v, t)    (v).x = (t); (v).y = 0.0
#define MAGMA_tally3_C_OP_NEG_ASGN(t, z) (t).x = -(z).x; (t).y = -(z).y
#define MAGMA_tally3_C_EQUAL(u,v)        (((u).x == (v).x) && ((u).y == (v).y))
#define MAGMA_tally3_C_GET_X(u)          ((u).x)
#define MAGMA_tally3_C_ASSIGN(v, t)      (v).x = (t).x; (v).y = (t).y
#define MAGMA_tally3_C_CNJG(v, t)        (v).x= (t).x; (v).y = -(t).y
#define MAGMA_tally3_C_DSCALE(v, t, s)   (v).x = (t).x/(s); (v).y = (t).y/(s)
#define MAGMA_tally3_C_OP_NEG(a, b, c)   (a).x = (b).x-(c).x; (a).y = (b).y-(c).y
#define MAGMA_tally3_C_MAKE(r, i)        make_cuFloatComplex((r), (i))
#define MAGMA_tally3_C_REAL(a)           cuCrealf(a)
#define MAGMA_tally3_C_IMAG(a)           cuCimagf(a)
#define MAGMA_tally3_C_ADD(a, b)         cuCaddf((a), (b))
#define MAGMA_tally3_C_SUB(a, b)         cuCsubf((a), (b))
#define MAGMA_tally3_C_MUL(a, b)         cuCmulf((a), (b))
#define MAGMA_tally3_C_DIV(a, b)         cuCdivf((a), (b))
#define MAGMA_tally3_C_ABS(a)            cuCabsf((a))
#define MAGMA_tally3_C_ZERO              make_cuFloatComplex(0.0, 0.0)
#define MAGMA_tally3_C_ONE               make_cuFloatComplex(1.0, 0.0)
#define MAGMA_tally3_C_HALF              make_cuFloatComplex(0.5, 0.0)
#define MAGMA_tally3_C_NEG_ONE           make_cuFloatComplex(-1.0, 0.0)
#define MAGMA_tally3_C_NEG_HALF          make_cuFloatComplex(-0.5, 0.0)

#define MAGMA_tally3_D_SET2REAL(v, t)    (v) = (t);
#define MAGMA_tally3_D_OP_NEG_ASGN(t, z) (t) = -(z)
#define MAGMA_tally3_D_EQUAL(u,v)        ((u) == (v))
#define MAGMA_tally3_D_GET_X(u)          (u)
#define MAGMA_tally3_D_ASSIGN(v, t)      (v) = (t)
#define MAGMA_tally3_D_CNJG(v, t)        (v) = (t)
#define MAGMA_tally3_D_DSCALE(v, t, s)   (v) = (t)/(s)
#define MAGMA_tally3_D_OP_NEG(a, b, c)   (a) = (b) - (c)
#define MAGMA_tally3_D_MAKE(r, i)        (r)
#define MAGMA_tally3_D_REAL(a)           (a)
#define MAGMA_tally3_D_IMAG(a)           (a)
#define MAGMA_tally3_D_ADD(a, b)         ( (a) + (b) )
#define MAGMA_tally3_D_SUB(a, b)         ( (a) - (b) )
#define MAGMA_tally3_D_MUL(a, b)         ( (a) * (b) )
#define MAGMA_tally3_D_DIV(a, b)         ( (a) / (b) )
#define MAGMA_tally3_D_ABS(a)            ((a)>0?(a):-(a))
#define MAGMA_tally3_D_ZERO              (0.0)
#define MAGMA_tally3_D_ONE               (1.0)
#define MAGMA_tally3_D_HALF              (0.5)
#define MAGMA_tally3_D_NEG_ONE           (-1.0)
#define MAGMA_tally3_D_NEG_HALF          (-0.5)

#define MAGMA_tally3_S_SET2REAL(v, t)    (v) = (t);
#define MAGMA_tally3_S_OP_NEG_ASGN(t, z) (t) = -(z)
#define MAGMA_tally3_S_EQUAL(u,v)        ((u) == (v))
#define MAGMA_tally3_S_GET_X(u)          (u)
#define MAGMA_tally3_S_ASSIGN(v, t)      (v) = (t)
#define MAGMA_tally3_S_CNJG(v, t)        (v) = (t)
#define MAGMA_tally3_S_DSCALE(v, t, s)   (v) = (t)/(s)
#define MAGMA_tally3_S_OP_NEG(a, b, c)   (a) = (b) - (c)
#define MAGMA_tally3_S_MAKE(r, i)        (r)
#define MAGMA_tally3_S_REAL(a)           (a)
#define MAGMA_tally3_S_IMAG(a)           (a)
#define MAGMA_tally3_S_ADD(a, b)         ( (a) + (b) )
#define MAGMA_tally3_S_SUB(a, b)         ( (a) - (b) )
#define MAGMA_tally3_S_MUL(a, b)         ( (a) * (b) )
#define MAGMA_tally3_S_DIV(a, b)         ( (a) / (b) )
#define MAGMA_tally3_S_ABS(a)            ((a)>0?(a):-(a))
#define MAGMA_tally3_S_ZERO              (0.0)
#define MAGMA_tally3_S_ONE               (1.0)
#define MAGMA_tally3_S_HALF              (0.5)
#define MAGMA_tally3_S_NEG_ONE           (-1.0)
#define MAGMA_tally3_S_NEG_HALF          (-0.5)

#ifndef CBLAS_SADDR
#define CBLAS_SADDR(a)  &(a)
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------
 *   -- MAGMA_tally3 function definitions
 * --------------------------------------------------------- */
void magma_tally3_xerbla( const char *name, magma_tally3_int_t info );
magma_tally3_context *magma_tally3_init(void *, void* (*func)(void *a), magma_tally3_int_t nthread, magma_tally3_int_t ncpu, 
                          magma_tally3_int_t ngpu, magma_tally3_int_t argc, char **argv);
void magma_tally3_finalize(magma_tally3_context *cntxt);
void auto_tune(char algorithm, char precision, magma_tally3_int_t ncores, magma_tally3_int_t ncorespsocket,
               magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t *nb, magma_tally3_int_t *ob, magma_tally3_int_t *ib,
               magma_tally3_int_t *nthreads, magma_tally3_int_t *nquarkthreads);



#ifdef __cplusplus
}
#endif

#endif

