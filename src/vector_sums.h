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

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/
/*---Struct declaration---*/

typedef struct {
  void* __restrict__ data;
  void* __restrict__ data_tmp; } GMVectorSums;

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

/*===========================================================================*/
/*---Compute---*/

/*---TODO: use these two functiuon only internally to the (pseudo-)class---*/

void gm_compute_float_vector_sums(GMVectors* vectors,
                                  GMFloat* __restrict__ vector_sums,
                                  GMFloat* __restrict__ vector_sums_tmp,
                                  GMEnv* env);

void gm_compute_bits2_vector_sums(GMVectors* vectors,
                                  GMFloat* __restrict__ vector_sums,
                                  GMFloat* __restrict__ vector_sums_tmp,
                                  GMEnv* env);

void GMVectorSums_compute(GMVectorSums* vector_sums,
                          GMVectors* vectors,
                          GMEnv* env);

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*===========================================================================*/

#endif /*---_vector_sums_h---*/

/*---------------------------------------------------------------------------*/
