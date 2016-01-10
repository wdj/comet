/**
 *
 * @file context.h
 *
 *  MAGMA_tally4 (version 1.6.1) --
 *  Univ. of Tennessee, Knoxville
 *  Univ. of California, Berkeley
 *  Univ. of Colorado, Denver
 *  @date January 2015
 *
 **/

#ifndef _MAGMA_tally4_CONTEXT_H_
#define _MAGMA_tally4_CONTEXT_H_


/* ------------------------------------------------------------
 * MAGMA_tally4 Context
 * --------------------------------------------------------- */
typedef struct magma_tally4_context_s
{
  /* Number of CPU core in this context */
  magma_tally4_int_t num_cores;

  /* Number of GPUs in this context */
  magma_tally4_int_t num_gpus;

  /* GPU contexts */
  CUcontext *gpu_context;

  /* QUARK scheduler */
  Quark *quark;

  /* Block size, internally used for some algorithms */
  magma_tally4_int_t nb;

  /* Pointer to other global algorithm-dependent parameters */ 
  void *params;

} magma_tally4_context_t;



magma_tally4_context_t * magma_tally4_init(void *, void* (*func)(void *a), magma_tally4_int_t nthread, magma_tally4_int_t ncpu, magma_tally4_int_t ngpu,
                magma_tally4_int_t argc, char **argv);
void magma_tally4_finalize(magma_tally4_context_t * cntxt);

void auto_tune(char algorithm, char precision, magma_tally4_int_t ncores, magma_tally4_int_t ncorespsocket,
               magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t *nb, magma_tally4_int_t *ob, magma_tally4_int_t *ib,
               magma_tally4_int_t *nthreads, magma_tally4_int_t *nquarkthreads);  


#endif /* _MAGMA_tally4_CONTEXT_H_ */