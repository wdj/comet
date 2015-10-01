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

#include "mpi.h"

#include "env.h"
#include "vectors.h"

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

  if (data_type_id == DATA_TYPE_ID_FLOAT) {
    vectors->data = malloc(num_field * num_vector_local * sizeof(Float_t));
  } else {
    Insist(Bool_false ? "Invalid data type." : 0);
  }
  vectors->data_type_id = data_type_id;
  vectors->num_field = num_field;
  vectors->num_vector_local = num_vector_local;

  int mpi_code;
  mpi_code = MPI_Allreduce(&(vectors->num_vector_local), &(vectors->num_vector),
                           1, MPI_INT, MPI_SUM, env->mpi_comm);
  Assert(mpi_code == MPI_SUCCESS);
  /*---FIX - check for int overflow---*/
  mpi_code = MPI_Allreduce(&(vectors->num_vector_local),
                           &(vectors->num_vector_local_max), 1, MPI_INT,
                           MPI_MAX, env->mpi_comm);
  Assert(mpi_code == MPI_SUCCESS);
}

/*===========================================================================*/
/*---Vectors pseudo-destructor---*/

void Vectors_destroy(Vectors* vectors, Env* env) {
  Assert(vectors);
  Assert(vectors->data);
  Assert(env);

  free(vectors->data);
  vectors->data = 0;
}

/*---------------------------------------------------------------------------*/
