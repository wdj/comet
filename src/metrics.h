/*---------------------------------------------------------------------------*/
/*!
 * \file   metrics.h
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Metrics pseudo-class, header.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

/*=============================================================================


=============================================================================*/

#ifndef _metrics_h_
#define _metrics_h_

#include "env.h"

/*===========================================================================*/
/*---Struct declaration---*/

typedef struct {
    int                num_vector;
    int                num_vector_local;
    int                num_vector_local_max;
    int                data_type_id;
    void* __restrict__ data;
} Metrics;

/*===========================================================================*/
/*---Functions---*/

void Metrics_create(Metrics* metrics, int data_type_id, int num_vector_local,
                    Env* env);

void Metrics_destroy(Metrics* metrics, Env* env);

/*---------------------------------------------------------------------------*/
/*---Accessors---*/





/*===========================================================================*/

#endif /*---_metrics_h_---*/

/*---------------------------------------------------------------------------*/
