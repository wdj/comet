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

void gm_compute_metrics(GMMetrics* metrics, GMVectors* vectors, GMEnv* env);

/*
GMFloat vector_sum(int len, GMFloat * const __restrict__ v1);
*/

/*===========================================================================*/

#endif /*---_compute_metrics_h---*/

/*---------------------------------------------------------------------------*/
