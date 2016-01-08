/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_metrics_2way.h
 * \author Wayne Joubert
 * \date   Thu Jan  7 10:21:09 EST 2016
 * \brief  Functions for computing 2-way metrics, header.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _compute_metrics_2way_h_
#define _compute_metrics_2way_h_

#include "env.h"
#include "vectors.h"
#include "metrics.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/

void gm_compute_metrics_2way_local(GMMetrics* metrics,
                                   GMVectors* vectors,
                                   GMEnv* env);

void gm_compute_metrics_2way_all2all(GMMetrics* metrics,
                                     GMVectors* vectors,
                                     GMEnv* env);

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*===========================================================================*/

#endif /*---_compute_metrics_2way_h_---*/

/*---------------------------------------------------------------------------*/
