/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_metrics_sorenson.c
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Functions for computing Sorenson metrics.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdlib.h>

#include "env.h"
#include "vectors.h"
#include "metrics.h"
#include "compute_metrics_sorenson.h"

#if 0
/*===========================================================================*/

void compute_metrics_sorenson_2way_cpu(Metrics* metrics,
                                       Vectors* vectors,
                                       Env* env) {
  Assert(metrics != NULL);
  Assert(vectors != NULL);
  Assert(env != NULL);

  Insist(env, ( ! env->all2all ) ? "Unimplemented." : 0);

  Insist(env, Bool_false ? "Unimplemented." : 0);
}

/*===========================================================================*/

void compute_metrics_sorenson_2way_gpu(Metrics* metrics,
                                       Vectors* vectors,
                                       Env* env) {
  Assert(metrics != NULL);
  Assert(vectors != NULL);
  Assert(env != NULL);

  Insist(env, ( ! env->all2all ) ? "Unimplemented." : 0);

  Insist(env, Bool_false ? "Unimplemented." : 0);
}

/*===========================================================================*/
#endif

/*---------------------------------------------------------------------------*/
