//-----------------------------------------------------------------------------
/*!
 * \file   vector_sums.hh
 * \author Wayne Joubert
 * \date   Sat Oct 24 10:48:32 EDT 2015
 * \brief  Per-vector computed quantities, header.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _gm_vector_sums_hh_
#define _gm_vector_sums_hh_

#include "env.hh"
#include "vectors.hh"

//=============================================================================
/*---Struct declaration---*/

typedef struct {
  GMFloat* __restrict__ sums;
  GMFloat* __restrict__ counts;
  GMFloat* __restrict__ sums_tmp_;
  GMFloat* __restrict__ counts_tmp_;
  size_t size_;
  size_t num_field_;
} GMVectorSums;

//=============================================================================
/*---Null object---*/

GMVectorSums GMVectorSums_null(void);

//=============================================================================
/*---Pseudo-constructor---*/

void GMVectorSums_create(GMVectorSums* this_,
                         int num_vector_local,
                         GMEnv* env);

//=============================================================================
/*---Pseudo-destructor---*/

void GMVectorSums_destroy(GMVectorSums* this_, GMEnv* env);

//=============================================================================
/*---Compute---*/

void GMVectorSums_compute(GMVectorSums* this_, GMVectors* vectors, GMEnv* env);

//=============================================================================
/*---Accessors---*/

GMFloat GMVectorSums_sum(const GMVectorSums* this_, int i,  GMEnv* env);

GMFloat GMVectorSums_count(const GMVectorSums* this_, int i,  GMEnv* env);

//-----------------------------------------------------------------------------

#endif // _gm_vector_sums_hh_

//-----------------------------------------------------------------------------
