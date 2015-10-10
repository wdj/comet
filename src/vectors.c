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
  memset((void*)&result, 0, sizeof(Vectors));
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
    case DATA_TYPE_ID_FLOAT: {
        vectors->num_field_dataval = num_field;
        vectors->data = malloc( vectors->num_field_dataval *
                                num_vector_local * sizeof(Float_t));
        Assert(vectors->data != NULL);
      } break;
    case DATA_TYPE_ID_BIT: {
        Assert( sizeof(Bits_t) == 8 );
        vectors->num_field_dataval = ceil_i( num_field, 8 * sizeof(Bits_t) );
        vectors->data = malloc( vectors->num_field_dataval *
                                num_vector_local * sizeof(Bits_t));
        Assert(vectors->data != NULL);
        /*---Ensure final pad bits of each vector are set to zero---*/
        int vector_local;
        for ( vector_local = 0; vector_local < num_vector_local;
              ++vector_local ) {
          const int dataval_num = Vectors_bit_dataval_num(vectors,
                            vectors->num_field_dataval-1, vector_local, env);
          ((ULInt_t*)vectors->data)[dataval_num] = 0;
        }
      } break;
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
  if (mpi_code) {
  } /*---Avoid unused variable warning---*/
  Assert(mpi_code == MPI_SUCCESS);

  size_t num_vector_bound = env->num_proc * (size_t)vectors->num_vector;
  if (num_vector_bound) {
  } /*---Avoid unused variable warning---*/
  Assert(num_vector_bound == (size_t)(int)num_vector_bound
             ? "Vector count too large to store in 32-bit int."
             : 0);

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
