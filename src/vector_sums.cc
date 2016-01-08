/*---------------------------------------------------------------------------*/
/*!
 * \file   vector_sums.c
 * \author Wayne Joubert
 * \date   Sat Oct 24 10:48:32 EDT 2015
 * \brief  Per-vector computed quantities.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdlib.h>

#include "env.h"
#include "vector_sums.h"
#include "vectors.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/
/*---Null object---*/

GMVectorSums GMVectorSums_null(void) {
  GMVectorSums x;
  x.data = NULL;
  x.data_tmp = NULL;
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

  switch (Env_metric_type(env)) {
    case GM_METRIC_TYPE_SORENSON: {
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    case GM_METRIC_TYPE_CZEKANOWSKI: {
      vector_sums->data = GMFloat_malloc(vectors->num_vector_local);
      vector_sums->data_tmp = Env_num_proc_field(env) == 1 ? NULL
                                   : GMFloat_malloc(vectors->num_vector_local);
    } break;
    case GM_METRIC_TYPE_CCC: {
      vector_sums->data = GMFloat_malloc(vectors->num_vector_local);
      vector_sums->data_tmp = Env_num_proc_field(env) == 1 ? NULL
                                   : GMFloat_malloc(vectors->num_vector_local);
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

  switch (Env_metric_type(env)) {
    case GM_METRIC_TYPE_SORENSON: {
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    case GM_METRIC_TYPE_CZEKANOWSKI: {
      GMAssert(vector_sums->data != NULL);
      free(vector_sums->data);
      vector_sums->data = NULL;
      if (vector_sums->data_tmp != NULL) {
        free(vector_sums->data_tmp);
        vector_sums->data_tmp = NULL;
      }
    } break;
    case GM_METRIC_TYPE_CCC: {
      GMAssert(vector_sums->data != NULL);
      free(vector_sums->data);
      vector_sums->data = NULL;
      if (vector_sums->data_tmp != NULL) {
        free(vector_sums->data_tmp);
        vector_sums->data_tmp = NULL;
      }
    } break;
    default:
      /*---Should never get here---*/
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
  } /*---case---*/
}

/*===========================================================================*/
/*---Compute the sum of elements of each vector on CPU, for denom---*/

void gm_compute_float_vector_sums(GMVectors* vectors,
                                  GMFloat* __restrict__ vector_sums,
                                  GMFloat* __restrict__ vector_sums_tmp,
                                  GMEnv* env) {
  GMAssert(vector_sums != NULL);
  GMAssert(vectors != NULL);
  GMAssert(vector_sums_tmp != NULL || Env_num_proc_field(env) == 1);
  GMAssert(env != NULL);

  const int num_proc = Env_num_proc_field(env);
  GMFloat* __restrict__ vector_sums_local = num_proc==1
                                          ? vector_sums : vector_sums_tmp;

  int i = 0;
  for (i = 0; i < vectors->num_vector_local; ++i) {
    GMFloat sum = 0;
    int field_local = 0;
    for (field_local = 0; field_local < vectors->num_field_local;
         ++field_local) {
      const GMFloat value = GMVectors_float_get(vectors, field_local, i, env);
      sum += value;
    }
    vector_sums_local[i] = sum;
  }

  /*---Do reduction across field procs if needed---*/

  if (num_proc > 1) {
    int mpi_code = 0;
    mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/
    mpi_code = MPI_Allreduce(vector_sums_local, vector_sums,
                 vectors->num_vector_local, GM_MPI_FLOAT, MPI_SUM,
                 Env_mpi_comm_field(env));
    GMAssert(mpi_code == MPI_SUCCESS);
  }
}

/*---------------------------------------------------------------------------*/

void gm_compute_bits2_vector_sums(GMVectors* vectors,
                                  GMFloat* __restrict__ vector_sums,
                                  GMFloat* __restrict__ vector_sums_tmp,
                                  GMEnv* env) {
  GMAssert(vector_sums != NULL);
  GMAssert(vectors != NULL);
  GMAssert(vector_sums_tmp != NULL || Env_num_proc_field(env) == 1);
  GMAssert(env != NULL);

  const int num_proc = Env_num_proc_field(env);
  GMFloat* __restrict__ vector_sums_local = num_proc==1
                                          ? vector_sums : vector_sums_tmp;

  /*----------*/
  if (env->compute_method_ == GM_COMPUTE_METHOD_REF) {
  /*----------*/
    int i = 0;
    for (i = 0; i < vectors->num_vector_local; ++i) {
      GMFloat sum = 0;
      int f = 0;
      for (f = 0; f < vectors->num_field_local; ++f) {
        /*---Slow way: sum each semi-nibble individually---*/
        const GMBits2 value = GMVectors_bits2_get(vectors, f, i, env);
        sum += ((value & 1) != 0) + ((value & 2) != 0);
      }
      vector_sums_local[i] = sum;
    }
  /*----------*/
  } else {
  /*----------*/
    int i = 0;
    for (i = 0; i < vectors->num_vector_local; ++i) {
      GMFloat sum = 0;
      int f = 0;
      for (f = 0; f < vectors->num_packedval_field_local; ++f) {
        /*---Fast way: sum all 64 bits of each word immediately---*/
        const GMBits2x64 value = GMVectors_bits2x64_get(vectors, f, i, env);
        sum += (GMFloat)gm_popcount64(value.data[0]);
        sum += (GMFloat)gm_popcount64(value.data[1]);
      }
      vector_sums_local[i] = sum;
    }
  /*----------*/
  } /*---if---*/
  /*----------*/

  /*---Do reduction across field procs if needed---*/

  if (num_proc > 1) {
    int mpi_code = 0;
    mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/
    mpi_code = MPI_Allreduce(vector_sums_local, vector_sums,
                 vectors->num_vector_local, GM_MPI_FLOAT, MPI_SUM,
                 Env_mpi_comm_field(env));
    GMAssert(mpi_code == MPI_SUCCESS);
  }
}

/*---------------------------------------------------------------------------*/

void GMVectorSums_compute(GMVectorSums* vector_sums,
                          GMVectors* vectors,
                          GMEnv* env) {
  GMAssert(vector_sums != NULL);
  GMAssert(vectors != NULL);
  GMAssert(env != NULL);

  switch (Env_metric_type(env)) {
    /*----------------------------------------*/
    case GM_METRIC_TYPE_SORENSON: {
      /*----------------------------------------*/

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    /*----------------------------------------*/
    case GM_METRIC_TYPE_CZEKANOWSKI: {
      /*----------------------------------------*/
      gm_compute_float_vector_sums(vectors, (GMFloat*)vector_sums->data,
                                        (GMFloat*)vector_sums->data_tmp,  env);
    } break;
    /*----------------------------------------*/
    case GM_METRIC_TYPE_CCC: {
      /*----------------------------------------*/
      gm_compute_bits2_vector_sums(vectors, (GMFloat*)vector_sums->data,
                                        (GMFloat*)vector_sums->data_tmp,  env);
    } break;
    /*----------------------------------------*/
    default:
      /*----------------------------------------*/
      /*---Should never get here---*/
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
  } /*---case---*/
}

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
