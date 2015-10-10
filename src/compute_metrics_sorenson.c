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

void gm_compute_metrics_sorenson_2way_cpu(GMMetrics* metrics,
                                       GMVectors* vectors,
                                       GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  GMInsist(env, ( ! env->all2all ) ? "Unimplemented." : 0);

  GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
}

/*===========================================================================*/

void gm_compute_metrics_sorenson_2way_gpu(GMMetrics* metrics,
                                       GMVectors* vectors,
                                       GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  GMInsist(env, ( ! env->all2all ) ? "Unimplemented." : 0);

  GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
}

/*===========================================================================*/
#endif

/*---------------------------------------------------------------------------*/
