/*---------------------------------------------------------------------------*/
/*!
 * \file   vector_sums.h
 * \author Wayne Joubert
 * \date   Sat Oct 24 10:48:32 EDT 2015
 * \brief  Per-vector computed quantities, header.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _vector_sums_h_
#define _vector_sums_h_

#include "env.h"
#include "vectors.h"

/*===========================================================================*/
/*---Struct declaration---*/

typedef struct {
  void* __restrict__ data;
} GMVectorSums;

/*===========================================================================*/
/*---Null object---*/

GMVectorSums GMVectorSums_null(void);

/*===========================================================================*/
/*---Pseudo-constructor---*/

void GMVectorSums_create(GMVectorSums* vector_sums,
                         GMVectors* vectors,
                         GMEnv* env);

/*===========================================================================*/
/*---Pseudo-destructor---*/

void GMVectorSums_destroy(GMVectorSums* vector_sums, GMEnv* env);

#endif /*---_vector_sums_h---*/

/*---------------------------------------------------------------------------*/

