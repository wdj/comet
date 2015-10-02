/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_metrics.h
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Functions for computing metrics, header.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _compute_metrics_h_
#define _compute_metrics_h_

#include "env.h"
#include "vectors.h"
#include "metrics.h"

/*===========================================================================*/

void compute_metrics(Metrics* metrics, Vectors* vectors, Env* env);
void compute_metrics_czek_2way_cpu(Metrics* metrics,
                                   Vectors* vectors,
                                   Env* env);
void compute_metrics_czek_2way_gpu(Metrics* metrics,
                                   Vectors* vectors,
                                   Env* env);

Float_t vector_sum(int len, Float_t * const __restrict__ v1);


/*===========================================================================*/

#endif /*---_compute_metrics_h---*/

/*---------------------------------------------------------------------------*/
