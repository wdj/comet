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

#include "mpi.h"

#include "env.h"
#include "metrics.h"

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

  if (env->global_all2all) {
    if (env->num_way == 2) {
      Insist(Bool_false ? "Not yet implemented." : 0); /*FIX*/
    } else /* (env->num_way == 3) */ {
      Insist(Bool_false ? "Not yet implemented." : 0); /*FIX*/
    }
  } else {
    if (env->num_way == 2) {
      /*---n choose 2---*/
      metrics->num_elts_local = (num_vector_local * (num_vector_local - 1)) / 2;
      metrics->index_map = malloc(metrics->num_elts_local * sizeof(int));
      int index = 0;
      int i = 0;
      for (i = 0; i < num_vector_local; ++i) {
        int j = 0;
        for (j = 0; j < i; ++j) {
          metrics->index_map[index++] = j + num_vector_local * i;
        }
      }
    } else /* (env->num_way == 3) */ {
      /*---n choose 3---*/
      metrics->num_elts_local =
          (num_vector_local * (num_vector_local - 1) * (num_vector_local - 2)) /
          6;
      metrics->index_map = malloc(metrics->num_elts_local * sizeof(int));
      int index = 0;
      int i = 0;
      for (i = 0; i < num_vector_local; ++i) {
        int j = 0;
        for (j = 0; j < i; ++j) {
          int k = 0;
          for (k = 0; k < j; ++k) {
            metrics->index_map[index++] =
                k + num_vector_local * (j + num_vector_local * i);
          }
        }
      }
    }
  }

  if (data_type_id == DATA_TYPE_ID_FLOAT) {
    metrics->data = malloc(metrics->num_elts_local * sizeof(Float_t));
    else if (data_type_id == DATA_TYPE_ID_BIT) {
      Insist(Bool_false ? "Not yet implemented." : 0); /*FIX*/
    }
    else {
      Insist(Bool_false ? "Invalid data type." : 0);
    }

    int mpi_code;
    mpi_code =
        MPI_Allreduce(&(metrics->num_vector_local), &(metrics->num_vector), 1,
                      MPI_INT, MPI_SUM, env->mpi_comm);
    Assert(mpi_code == MPI_SUCCESS);
    mpi_code = MPI_Allreduce(&(metrics->num_vector_local),
                             &(metrics->num_vector_local_max), 1, MPI_INT,
                             MPI_MAX, env->mpi_comm);
    Assert(mpi_code == MPI_SUCCESS);
  }

  /*===========================================================================*/
  /*---Metrics pseudo-destructor---*/

  void Metrics_destroy(Metrics * metrics, Env * env) {
    Assert(metrics);
    Assert(metrics->data);
    Assert(env);

    free(metrics->data);
    free(metrics->index_map);
    metrics->data = 0;
  }

/*---------------------------------------------------------------------------*/
