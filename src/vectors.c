/*---------------------------------------------------------------------------*/
/*!
 * \file   vectors.c
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Vectors pseudo-class.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>

#include "mpi.h"

#include "env.h"
#include "vectors.h"

/*===========================================================================*/
/*---Null object---*/

Vectors Vectors_null() {
  Vectors result;
  memset( (void*)&result, 0, sizeof(Vectors) );
  return result;
}

/*===========================================================================*/
/*---Vectors pseudo-constructor---*/

void Vectors_create(Vectors* vectors,
                    int data_type_id,
                    int num_field,
                    int num_vector_local,
                    Env* env) {
  Assert(vectors);
  Assert(num_field >= 0);
  Assert(num_vector_local >= 0);
  Assert(env);

  /*---Allocations---*/
  switch (data_type_id) {
    case DATA_TYPE_ID_FLOAT:
      vectors->data = malloc(num_field * num_vector_local * sizeof(Float_t));
      Assert(vectors->data != NULL);
      break;
    case DATA_TYPE_ID_BIT:
      Insist(env, Bool_false ? "Unimplemented." : 0);
      break;
    default:
      Assert(Bool_false ? "Invalid data type." : 0);
  }

  vectors->data_type_id = data_type_id;
  vectors->num_field = num_field;
  vectors->num_vector_local = num_vector_local;

  /*---Compute global values---*/

  int mpi_code = MPI_Allreduce(&(vectors->num_vector_local),
                               &(vectors->num_vector_local_max), 1, MPI_INT,
                               MPI_MAX, env->mpi_comm);
  Assert(mpi_code == MPI_SUCCESS);

  size_t num_vector_bound = env->num_proc * (size_t) vectors->num_vector;
  Assert( num_vector_bound == (size_t)(int)num_vector_bound
            ? "Vector count too large to store in 32-bit int." : 0 );

  mpi_code = MPI_Allreduce(&(vectors->num_vector_local), &(vectors->num_vector),
                           1, MPI_INT, MPI_SUM, env->mpi_comm);
  Assert(mpi_code == MPI_SUCCESS);
}

/*===========================================================================*/
/*---Vectors pseudo-destructor---*/

void Vectors_destroy(Vectors* vectors, Env* env) {
  Assert(vectors);
  Assert(vectors->data);
  Assert(env);

  free(vectors->data);
  *vectors = Vectors_null();
}

/*---------------------------------------------------------------------------*/
