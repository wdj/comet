/*---------------------------------------------------------------------------*/
/*!
 * \file   metrics.c
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Metrics pseudo-class.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>

#include "mpi.h"

#include "env.h"
#include "metrics.h"

/*===========================================================================*/
/*---Null object---*/

GMMetrics GMMetrics_null() {
  GMMetrics result;
  memset((void*)&result, 0, sizeof(GMMetrics));
  return result;
}

/*===========================================================================*/
/*---Metrics pseudo-constructor---*/

void GMMetrics_create(GMMetrics* metrics,
                    int data_type_id,
                    int num_vector_local,
                    GMEnv* env) {
  GMAssert(metrics);
  GMAssert(num_vector_local >= 0);
  GMAssert(env);

  metrics->data_type_id = data_type_id;
  metrics->num_vector_local = num_vector_local;

  /*---Compute global values---*/

  int mpi_code = MPI_Allreduce(&(metrics->num_vector_local),
                               &(metrics->num_vector_local_max), 1, MPI_INT,
                               MPI_MAX, env->mpi_comm);
  if (mpi_code) {
  } /*---Avoid unused variable warning---*/
  GMAssert(mpi_code == MPI_SUCCESS);

  size_t num_vector_bound = env->num_proc * (size_t)metrics->num_vector;
  if (num_vector_bound) {
  } /*---Avoid unused variable warning---*/
  GMAssert(num_vector_bound == (size_t)(int)num_vector_bound
             ? "Vector count too large to store in 32-bit int."
             : 0);

  mpi_code = MPI_Allreduce(&(metrics->num_vector_local), &(metrics->num_vector),
                           1, MPI_INT, MPI_SUM, env->mpi_comm);
  GMAssert(mpi_code == MPI_SUCCESS);

  /*---Compute number of elements etc.---*/

  if (env->all2all) {
    if (env->num_way == 2) {
      /*---Store strict upper triang of diag block and half
          ` the off-diag blocks---*/
      metrics->num_elts_local = gm_nchoosek(num_vector_local, env->num_way) +
          ( env->num_proc / 2 ) * num_vector_local * num_vector_local;
      metrics->coords_global_from_index = malloc( metrics->num_elts_local *
                                                  sizeof(size_t));
      GMAssert(metrics->coords_global_from_index != NULL);
      int index = 0;
      int i = 0;

/*FIX*/
      for (i = 0; i < num_vector_local; ++i) {
        const size_t i_global = i + num_vector_local * env->proc_num;
        /*---j here is a global index---*/
        const int beg = env->proc_num * num_vector_local + i + 1;
        const int end = ( env->proc_num + 1 + ( env->num_proc / 2 ) )
                        * num_vector_local;
        size_t j_global_unwrapped = 0;
        for (j_global_unwrapped = beg; j_global_unwrapped < end;
             ++j_global_unwrapped) {
          const size_t j_global = j_global_unwrapped % metrics->num_vector;
          metrics->coords_global_from_index[index++] = i_global +
                                               metrics->num_vector * j_global;
        }
      } /*---i---*/
    } else /* (env->num_way == 3) */ {

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    }
  } else { /*---if not all2all---*/
    metrics->num_elts_local = gm_nchoosek(num_vector_local, env->num_way);
    metrics->coords_global_from_index = malloc( metrics->num_elts_local *
                                                sizeof(size_t));
    GMAssert(metrics->coords_global_from_index != NULL);
    /*---LATER: generalize this to N-way---*/
    if (env->num_way == 2) {
      /*---Need store only strict upper triangular part of matrix---*/
      int index = 0;
      int j = 0;
      for (j = 0; j < num_vector_local; ++j) {
        const size_t j_global = j + num_vector_local * env->proc_num;
        int i = 0;
        for (i = 0; i < j; ++i) {
          const size_t i_global = i + num_vector_local * env->proc_num;
          metrics->coords_global_from_index[index++] =
              i_global + metrics->num_vector * j_global;
        }
      }
    } else /* (env->num_way == 3) */ {
      /*---Need store only strict interior of tetrahedron---*/
      int index = 0;
      int k = 0;
      for (k = 0; k < num_vector_local; ++k) {
        const size_t k_global = k + num_vector_local * env->proc_num;
        int j = 0;
        for (j = 0; j < k; ++j) {
          const size_t j_global = j + num_vector_local * env->proc_num;
          int i = 0;
          for (i = 0; i < j; ++k) {
            const size_t i_global = i + num_vector_local * env->proc_num;
            metrics->coords_global_from_index[index++] =
                i_global + metrics->num_vector * (
                j_global + metrics->num_vector * k_global );
          }
        }
      }
    } /*---if num_way---*/
  }   /*---if all2all---*/

  /*---Allocations---*/

  switch (data_type_id) {
    case GM_DATA_TYPE_FLOAT:
      metrics->data = malloc(metrics->num_elts_local * sizeof(GMFloat));
      GMAssert(metrics->data != NULL);
      break;
    case GM_DATA_TYPE_BIT: {
      const size_t num_floats_needed = gm_ceil_i8( metrics->num_elts_local,
                                                   8 * sizeof(GMFloat) );
      metrics->data = malloc(num_floats_needed * sizeof(GMFloat));
      GMAssert(metrics->data != NULL);
      } break;
    default:
      GMAssert(GM_BOOL_FALSE ? "Invalid data type." : 0);
  }
}

/*===========================================================================*/
/*---Metrics pseudo-destructor---*/

void GMMetrics_destroy(GMMetrics* metrics, GMEnv* env) {
  GMAssert(metrics);
  GMAssert(metrics->data);
  GMAssert(env);

  free(metrics->data);
  free(metrics->coords_global_from_index);
  *metrics = GMMetrics_null();
}

/*===========================================================================*/
/*---Metrics checksum---*/

/* This should be invariant, up to roundoff, on CPU vs. GPU. */

double GMMetrics_checksum(GMMetrics* metrics, GMEnv* env) {
  GMAssert(metrics);
  GMAssert(metrics->data);
  GMAssert(env);

  double result = 0;

  if (env->all2all) {
    GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
  } else {
    switch (metrics->data_type_id) {
      case GM_DATA_TYPE_FLOAT: {
        int i = 0;
        for (i = 0; i < metrics->num_elts_local; ++i) {
          const size_t i_global =
              i + metrics->num_elts_local * (size_t)env->num_proc;
          result += ((GMFloat*)metrics->data)[i] * gm_randomize(i_global);
        } /*---for i---*/
      } break;
      case GM_DATA_TYPE_BIT:
        GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
        break;
      default:
        GMAssert(GM_BOOL_FALSE ? "Invalid data type." : 0);
    } /*---switch---*/
    const double tmp = result;
    int mpi_code =
        MPI_Allreduce(&tmp, &result, 1, MPI_DOUBLE, MPI_MAX, env->mpi_comm);
    if (mpi_code) {
    } /*---Avoid unused variable warning---*/
    GMAssert(mpi_code == MPI_SUCCESS);
  } /*---if all2all---*/

  return result;
}

/*---------------------------------------------------------------------------*/
