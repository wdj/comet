/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_metrics_czekanowski_2way.h
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Functions for computing 2-way Czekanowski metrics, header.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _compute_metrics_czekanowski_2way_h_
#define _compute_metrics_czekanowski_2way_h_

#include "env.h"
#include "vectors.h"
#include "metrics.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/

void gm_compute_metrics_czekanowski_2way_cpu(GMMetrics* metrics,
                                             GMVectors* vectors,
                                             GMEnv* env);

void gm_compute_metrics_czekanowski_2way_gpu(GMMetrics* metrics,
                                             GMVectors* vectors,
                                             GMEnv* env);

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*===========================================================================*/

#endif /*---_compute_metrics_czekanowski_2way_h_---*/

/*---------------------------------------------------------------------------*/
