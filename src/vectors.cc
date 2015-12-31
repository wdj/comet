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

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/
/*---Null object---*/

GMVectors GMVectors_null() {
  GMVectors result;
  memset((void*)&result, 0, sizeof(GMVectors));
  return result;
}

/*===========================================================================*/
/*---Vectors pseudo-constructor---*/

void GMVectors_create(GMVectors* vectors,
                      int data_type_id,
                      int num_field,
                      int num_vector_local,
                      GMEnv* env) {
  GMAssert(vectors);
  GMAssert(num_field >= 0);
  GMAssert(num_vector_local >= 0);
  GMAssert(env);

  *vectors = GMVectors_null();

  if (!Env_is_proc_active(env)) {
    return;
  }

  GMInsist(env, num_field % Env_num_proc_field(env) == 0
    ? "num_proc_field must exactly divide the total number of fields" : 0);

  vectors->data_type_id = data_type_id;
  vectors->num_field = num_field;
  vectors->num_field_local = num_field / Env_num_proc_field(env);
  vectors->num_vector_local = num_vector_local;

  /*---Compute global values---*/

  size_t num_vector_bound = 0;
  num_vector_bound =
      num_vector_bound * 1; /*---Avoid unused variable warning---*/
  num_vector_bound = Env_num_proc_vector(env) * (size_t)num_vector_local;
  GMAssert(num_vector_bound == (size_t)(int)num_vector_bound
               ? "Vector count too large to store in 32-bit int."
               : 0);

  int mpi_code = 0;
  mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/
  mpi_code = MPI_Allreduce(&(vectors->num_vector_local), &(vectors->num_vector),
                           1, MPI_INT, MPI_SUM, Env_mpi_comm_vector(env));
  GMAssert(mpi_code == MPI_SUCCESS);

  /*---Set element sizes---*/

  switch (data_type_id) {
    /*--------------------*/
    case GM_DATA_TYPE_FLOAT: {
      vectors->num_bits_per_val = 8 * sizeof(GMFloat);
      vectors->num_bits_per_packedval = 8 * sizeof(GMFloat);
    } break;
    /*--------------------*/
    case GM_DATA_TYPE_BITS2: {
      vectors->num_bits_per_val = 2;
      vectors->num_bits_per_packedval = 8 * sizeof(GMBits2x64);
    } break;
    /*--------------------*/
    case GM_DATA_TYPE_BITS1: {
      //---(design is not complete)
      vectors->num_bits_per_val = 1;
      vectors->num_bits_per_packedval = 8 * sizeof(GMBits1x64);
    } break;
    /*--------------------*/
    default:
      GMAssert(GM_BOOL_FALSE ? "Invalid data type." : 0);
  }

  /*---Allocations---*/

  vectors->num_packedval_field_local =
    gm_ceil_i8(vectors->num_field_local *
               (size_t)(vectors->num_bits_per_val),
               vectors->num_bits_per_packedval);
  vectors->num_packedval_local =
      vectors->num_packedval_field_local * (size_t)num_vector_local;
  vectors->data = malloc(vectors->num_packedval_local *
                         (vectors->num_bits_per_packedval / 8 ));
  GMAssert(vectors->data != NULL);


  /*---Ensure final pad bits of each vector are set to zero---*/

  switch (data_type_id) {
    /*--------------------*/
    case GM_DATA_TYPE_FLOAT: {
      /*---NO-OP---*/
    } break;
    /*--------------------*/
    case GM_DATA_TYPE_BITS2: {
      int vector_local;
      if (vectors->num_field_local > 0) {
        const int packedval_field_local
                                     = vectors->num_packedval_field_local - 1;
        GMBits2x64 zero = GMBits2x64_null();
        for (vector_local=0; vector_local<num_vector_local; ++vector_local) {
          GMVectors_bits2x64_set(vectors, packedval_field_local, vector_local,
                                 zero, env);
        }
      }
    } break;
    /*--------------------*/
    case GM_DATA_TYPE_BITS1: {
      int vector_local;
      if (vectors->num_field_local > 0) {
        const int packedval_field_local
                                     = vectors->num_packedval_field_local - 1;
        GMBits1x64 zero = GMBits1x64_null();
        for (vector_local=0; vector_local<num_vector_local; ++vector_local) {
          GMVectors_bits1x64_set(vectors, packedval_field_local, vector_local,
                                 zero, env);
        }
      }
    } break;
    /*--------------------*/
    default:
      GMAssert(GM_BOOL_FALSE ? "Invalid data type." : 0);
  }

}

/*===========================================================================*/
/*---Vectors pseudo-destructor---*/

void GMVectors_destroy(GMVectors* vectors, GMEnv* env) {
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);
  GMAssert(vectors->data != NULL || !Env_is_proc_active(env));

  if (!Env_is_proc_active(env)) {
    return;
  }

  free(vectors->data);
  *vectors = GMVectors_null();
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
