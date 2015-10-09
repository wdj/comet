/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_metrics_utils.h
 * \author Wayne Joubert
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Functions for computing metrics, utilities, header.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _compute_metrics_utils_h_
#define _compute_metrics_utils_h_

#include "env.h"
#include "vectors.h"
#include "metrics.h"

/*===========================================================================*/

void compute_vector_sums(Vectors* vectors,
                         Float_t* __restrict__ vector_sums,
                         Env* env);

/*===========================================================================*/

#endif /*---_compute_metrics_utils_h---*/

/*---------------------------------------------------------------------------*/
