/*---------------------------------------------------------------------------*/
/*!
 * \file   vector_sums.c
 * \author Wayne Joubert
 * \date   Sat Oct 24 10:48:32 EDT 2015
 * \brief  Per-vector computed quantities.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include "env.h"
#include "vector_sums.h"
#include "vectors.h"

/*===========================================================================*/
/*---Null object---*/

GMVectorSums GMVectorSums_null(void) {
  GMVectorSums x;
  x.data = NULL;
  return x;
}

/*===========================================================================*/
/*---Pseudo-constructor---*/

void GMVectorSums_create(GMVectorSums* vector_sums,
                         GMVectors* vectors,
                         GMEnv* env) {
  GMAssert(vector_sums);
  GMAssert(vectors);
  GMAssert(env);

  switch (env->metric_type) {
    case GM_METRIC_TYPE_SORENSON: {

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    case GM_METRIC_TYPE_CZEKANOWSKI: {

      vector_sums->data = GMFloat_malloc(vectors->num_vector_local);

    } break;
    case GM_METRIC_TYPE_CCC: {

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    default:
      /*---Should never get here---*/
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
  } /*---case---*/

}

/*===========================================================================*/
/*---Pseudo-destructor---*/

void GMVectorSums_destroy(GMVectorSums* vector_sums, GMEnv* env) {
  GMAssert(vector_sums);
  GMAssert(env);

  switch (env->metric_type) {
    case GM_METRIC_TYPE_SORENSON: {

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    case GM_METRIC_TYPE_CZEKANOWSKI: {

      GMAssert(vector_sums->data != NULL);
      free(vector_sums->data);
      vector_sums->data = NULL;

    } break;
    case GM_METRIC_TYPE_CCC: {

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    default:
      /*---Should never get here---*/
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
  } /*---case---*/

}

/*---------------------------------------------------------------------------*/

