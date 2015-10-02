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

Metrics Metrics_null() {
  Metrics result;
  memset( (void*)&result, 0, sizeof(Metrics) );
  return result;
}

/*===========================================================================*/
/*---Metrics pseudo-constructor---*/

void Metrics_create(Metrics* metrics,
                    int data_type_id,
                    int num_vector_local,
                    Env* env) {
  Assert(metrics);
  Assert(num_vector_local >= 0);
  Assert(env);

  metrics->data_type_id = data_type_id;
  metrics->num_vector_local = num_vector_local;

  /*---Compute global values---*/

  int mpi_code = MPI_Allreduce(&(metrics->num_vector_local),
                               &(metrics->num_vector_local_max), 1, MPI_INT,
                               MPI_MAX, env->mpi_comm);
  Assert(mpi_code == MPI_SUCCESS);

  size_t num_vector_bound = env->num_proc * (size_t) metrics->num_vector;
  Assert( num_vector_bound == (size_t)(int)num_vector_bound
            ? "Vector count too large to store in 32-bit int." : 0 );

  mpi_code = MPI_Allreduce(&(metrics->num_vector_local), &(metrics->num_vector),
                           1, MPI_INT, MPI_SUM, env->mpi_comm);
  Assert(mpi_code == MPI_SUCCESS);

  /*---Compute number of elements etc.---*/

  if (env->global_all2all) {
    if (env->num_way == 2) {
      Insist(env, Bool_false ? "Unimplemented." : 0);
    } else /* (env->num_way == 3) */ {
      Insist(env, Bool_false ? "Unimplemented." : 0);
    }
  } else {
    metrics->num_elts_local = nchoosek( num_vector_local, env->num_way );
    metrics->index_map = malloc(metrics->num_elts_local * sizeof(size_t));
    Assert(metrics->index_map != NULL);
    /*---TODO: generalize this to N-way---*/
    if (env->num_way == 2) {
      int index = 0;
      int i = 0;
      for (i = 0; i < num_vector_local; ++i) {
        int j = 0;
        for (j = i+1; j < num_vector_local; ++j) {
          metrics->index_map[index++] = i + num_vector_local * (size_t)j;
        }
      }
    } else /* (env->num_way == 3) */ {
      int index = 0;
      int i = 0;
      for (i = 0; i < num_vector_local; ++i) {
        int j = 0;
        for (j = i+1; j < num_vector_local; ++j) {
          int k = 0;
          for (k = j+1; k < num_vector_local; ++k) {
            metrics->index_map[index++] =
                i + num_vector_local * (j + num_vector_local * (size_t)k);
          }
        }
      }
    } /*---if num_way---*/
  } /*---if global_all2all---*/

  /*---Allocations---*/

  switch (data_type_id) {
    case DATA_TYPE_ID_FLOAT:
      metrics->data = malloc(metrics->num_elts_local * sizeof(Float_t));
      Assert(metrics->data != NULL);
      break;
    case DATA_TYPE_ID_BIT:
      Insist(env, Bool_false ? "Unimplemented." : 0);
      break;
    default:
      Assert(Bool_false ? "Invalid data type." : 0);
  }
}

/*===========================================================================*/
/*---Metrics pseudo-destructor---*/

void Metrics_destroy(Metrics * metrics, Env * env) {
  Assert(metrics);
  Assert(metrics->data);
  Assert(env);

  free(metrics->data);
  free(metrics->index_map);
  *metrics = Metrics_null();
}

/*===========================================================================*/
/*---Metrics checksum---*/

/* This should be invariant, up to roundoff, on CPU vs. GPU. */

double Metrics_checksum(Metrics * metrics, Env * env) {
  Assert(metrics);
  Assert(metrics->data);
  Assert(env);

  double result = 0;

  if ( env->global_all2all ) {
    Insist(env, Bool_false ? "Unimplemented." : 0);
  } else {
    switch (metrics->data_type_id) {
      case DATA_TYPE_ID_FLOAT:
        {
          int i = 0;
          for ( i = 0; i < metrics->num_elts_local; ++i ) {
            const size_t i_global = i + metrics->num_elts_local
                                      * (size_t) env->num_proc;
            result += ((Float_t*)metrics->data)[i] * randomize( i_global );
          } /*---for i---*/
        }
        break;
      case DATA_TYPE_ID_BIT:
        Insist(env, Bool_false ? "Unimplemented." : 0);
        break;
      default:
        Assert(Bool_false ? "Invalid data type." : 0);
    } /*---switch---*/
    const double tmp = result;
    int mpi_code = MPI_Allreduce(&tmp, &result, 1, MPI_DOUBLE,
                                 MPI_MAX, env->mpi_comm);
    Assert(mpi_code == MPI_SUCCESS);
  } /*---if global_all2all---*/

  return result;
}

/*---------------------------------------------------------------------------*/
