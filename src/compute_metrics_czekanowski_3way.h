/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_metrics_czekanowski_3way.h
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Functions for computing 3-way Czekanowski metrics, header.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _compute_metrics_czekanowski_3way_h_
#define _compute_metrics_czekanowski_3way_h_

#include "env.h"
#include "vectors.h"
#include "metrics.h"

/*===========================================================================*/

void compute_metrics_czekanowski_3way_cpu(Metrics* metrics,
                                          Vectors* vectors,
                                          Env* env);
void compute_metrics_czekanowski_3way_gpu(Metrics* metrics,
                                          Vectors* vectors,
                                          Env* env);

/*===========================================================================*/

#endif /*---_compute_metrics_czekanowski_3way_h_---*/

/*---------------------------------------------------------------------------*/
