/**
 *
 * @file context.h
 *
 *  MAGMA_minproduct (version 1.6.1) --
 *  Univ. of Tennessee, Knoxville
 *  Univ. of California, Berkeley
 *  Univ. of Colorado, Denver
 *  @date January 2015
 *
 **/

#ifndef _MAGMA_minproduct_CONTEXT_H_
#define _MAGMA_minproduct_CONTEXT_H_


/* ------------------------------------------------------------
 * MAGMA_minproduct Context
 * --------------------------------------------------------- */
typedef struct magma_minproduct_context_s
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

} magma_minproduct_context_t;



magma_minproduct_context_t * magma_minproduct_init(void *, void* (*func)(void *a), magma_minproduct_int_t nthread, magma_minproduct_int_t ncpu, magma_minproduct_int_t ngpu,
                magma_minproduct_int_t argc, char **argv);
void magma_minproduct_finalize(magma_minproduct_context_t * cntxt);

void auto_tune(char algorithm, char precision, magma_minproduct_int_t ncores, magma_minproduct_int_t ncorespsocket,
               magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t *nb, magma_minproduct_int_t *ob, magma_minproduct_int_t *ib,
               magma_minproduct_int_t *nthreads, magma_minproduct_int_t *nquarkthreads);  


#endif /* _MAGMA_minproduct_CONTEXT_H_ */
