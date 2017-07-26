/*---------------------------------------------------------------------------*/
/*!
 * \file   vector_sums.cc
 * \author Wayne Joubert
 * \date   Sat Oct 24 10:48:32 EDT 2015
 * \brief  Per-vector computed quantities.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdlib.h>

#include "env.hh"
#include "vector_sums.hh"
#include "vectors.hh"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/
/*---Null object---*/

GMVectorSums GMVectorSums_null(void) {
  GMVectorSums x;
  x.sums = NULL;
  x.counts = NULL;
  x.sums_tmp = NULL;
  x.counts_tmp = NULL;
  x.size = 0;
  return x;
}

/*===========================================================================*/
/*---Pseudo-constructor---*/

void GMVectorSums_create(GMVectorSums* vector_sums,
                         GMVectors* vectors,
                         GMEnv* env) {
  GMAssertAlways(vector_sums);
  GMAssertAlways(vectors);
  GMAssertAlways(env);

  vector_sums->size = vectors->num_vector_local;
  const int num_proc = GMEnv_num_proc_field(env);

  switch (GMEnv_metric_type(env)) {
    case GM_METRIC_TYPE_SORENSON: {
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    case GM_METRIC_TYPE_CZEKANOWSKI: {
      vector_sums->sums = GMFloat_malloc(vectors->num_vector_local, env);
      vector_sums->sums_tmp = num_proc == 1
                              ? NULL
                              : GMFloat_malloc(vectors->num_vector_local, env);
    } break;
    case GM_METRIC_TYPE_CCC: {
      vector_sums->sums = GMFloat_malloc(vectors->num_vector_local, env);
      vector_sums->sums_tmp = num_proc == 1
                              ? NULL
                              : GMFloat_malloc(vectors->num_vector_local, env);
      if (env->sparse) {
        vector_sums->counts = GMFloat_malloc(vectors->num_vector_local, env);
        vector_sums->counts_tmp = num_proc == 1
                              ? NULL
                              : GMFloat_malloc(vectors->num_vector_local, env);
      }
    } break;
    default:
      /*---Should never get here---*/
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
  } /*---case---*/
}

/*===========================================================================*/
/*---Pseudo-destructor---*/

void GMVectorSums_destroy(GMVectorSums* vector_sums, GMEnv* env) {
  GMAssertAlways(vector_sums);
  GMAssertAlways(env);

  switch (GMEnv_metric_type(env)) {
    case GM_METRIC_TYPE_SORENSON: {
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    case GM_METRIC_TYPE_CZEKANOWSKI: {
      GMAssertAlways(vector_sums->sums != NULL);
      GMFloat_free((GMFloat*)vector_sums->sums, vector_sums->size, env);
      vector_sums->sums = NULL;
      if (vector_sums->sums_tmp != NULL) {
        GMFloat_free((GMFloat*)vector_sums->sums_tmp, vector_sums->size, env);
        vector_sums->sums_tmp = NULL;
      }
    } break;
    case GM_METRIC_TYPE_CCC: {
      GMAssertAlways(vector_sums->sums != NULL);
      GMFloat_free((GMFloat*)vector_sums->sums, vector_sums->size, env);
      vector_sums->sums = NULL;
      if (vector_sums->sums_tmp != NULL) {
        GMFloat_free((GMFloat*)vector_sums->sums_tmp, vector_sums->size, env);
        vector_sums->sums_tmp = NULL;
      }
      if (vector_sums->counts_tmp != NULL) {
        GMFloat_free((GMFloat*)vector_sums->counts_tmp, vector_sums->size, env);
        vector_sums->counts_tmp = NULL;
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
                                  GMFloat* __restrict__ sums,
                                  GMFloat* __restrict__ sums_tmp,
                                  GMEnv* env) {
  GMAssertAlways(vectors && sums && env);
  GMAssertAlways(sums_tmp != NULL || GMEnv_num_proc_field(env) == 1);

  const int num_proc = GMEnv_num_proc_field(env);
  GMFloat* const sums_local = num_proc == 1 ? sums : sums_tmp;

  /*---Sum up all values in each vector---*/

  int i = 0;
#pragma omp parallel for
  for (i = 0; i < vectors->num_vector_local; ++i) {
    GMFloat sum = 0;
    int fl = 0;
    //#pragma omp parallel for reduction(+:sum)
    for (fl = 0; fl < vectors->num_field_local; ++fl) {
      const GMFloat value = GMVectors_float_get(vectors, fl, i, env);
      sum += value;
    }
    sums_local[i] = sum;
  }

  /*---Do reduction across field procs if needed---*/

  if (num_proc > 1) {
    int mpi_code = 0;
    mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/
    mpi_code =
        MPI_Allreduce(sums_local, sums, vectors->num_vector_local,
                      GM_MPI_FLOAT, MPI_SUM, GMEnv_mpi_comm_field(env));
    GMAssertAlways(mpi_code == MPI_SUCCESS);
  }

  env->ops_local += 2 * vectors->num_vector_local *
                    (double)vectors->num_field_local;
}

/*---------------------------------------------------------------------------*/

void gm_compute_bits2_vector_sums(GMVectors* vectors,
                                  GMFloat* __restrict__ sums,
                                  GMFloat* __restrict__ sums_tmp,
                                  GMFloat* __restrict__ counts,
                                  GMFloat* __restrict__ counts_tmp,
                                  GMEnv* env) {
  GMAssertAlways(vectors && sums && env);
  GMAssertAlways(sums_tmp != NULL || GMEnv_num_proc_field(env) == 1);
  GMAssertAlways(counts_tmp != NULL || GMEnv_num_proc_field(env) == 1 ||
                 !env->sparse);

  const int num_proc = GMEnv_num_proc_field(env);
  GMFloat* const sums_local = num_proc == 1 ? sums : sums_tmp;
  GMFloat* const counts_local = num_proc == 1 ? counts : counts_tmp;

  /*---Count number of 1-bits in each vector---*/

  /*----------*/
  if (env->compute_method_ == GM_COMPUTE_METHOD_REF) {
    /*----------*/
    int i = 0;
#pragma omp parallel for
    for (i = 0; i < vectors->num_vector_local; ++i) {
      GMFloat sum = 0;
      int fl = 0;
      if (env->sparse) {
        GMFloat count = 0;
        for (fl = 0; fl < vectors->num_field_local; ++fl) {
          /*---Slow way: sum each semi-nibble individually---*/
          const GMBits2 value = GMVectors_bits2_get(vectors, fl, i, env);
          if (!( ((value & 1) == 0) && ((value & 2) != 0) )){
            sum += ((value & 1) != 0) + ((value & 2) != 0);
            count++;
          }
        }
        GMAssert(sum >= 0 && sum <= 2 * vectors->num_field);
        GMAssert(count >= 0 && count <= vectors->num_field);
        sums_local[i] = sum;
        counts_local[i] = count;
      } else { /*---sparse---*/
        //#pragma omp parallel for reduction(+:sum)
        for (fl = 0; fl < vectors->num_field_local; ++fl) {
          /*---Slow way: sum each semi-nibble individually---*/
          const GMBits2 value = GMVectors_bits2_get(vectors, fl, i, env);
          sum += ((value & 1) != 0) + ((value & 2) != 0);
        }
        GMAssert(sum >= 0 && sum <= 2 * vectors->num_field);
        sums_local[i] = sum;
      } /*---sparse---*/
    } /*---for i---*/
    /*----------*/
  } else {
    /*----------*/
    int i = 0;
    for (i = 0; i < vectors->num_vector_local; ++i) {
      GMFloat sum = 0;
      int f = 0;
      if (env->sparse) {
        const GMUInt64 oddbits = 0x5555555555555555;
        GMFloat count = 0;
        for (f = 0; f < vectors->num_packedval_field_local; ++f) {
          /*---Fast way: sum all 64 bits of each word immediately---*/
          const GMBits2x64 value = GMVectors_bits2x64_get(vectors, f, i, env);
          const GMUInt64 data0 = value.data[0];
          const GMUInt64 data1 = value.data[1];
          const GMUInt64 oddmask0 = (data0 | ~(data0 >> 1)) & oddbits;
          const GMUInt64 oddmask1 = (data1 | ~(data1 >> 1)) & oddbits;
          const GMUInt64 mask0 = oddmask0 | (oddmask0 << 1);
          const GMUInt64 mask1 = oddmask1 | (oddmask1 << 1);
          sum += (GMFloat)gm_popcount64(data0 & mask0);
          sum += (GMFloat)gm_popcount64(data1 & mask1);
          count += (GMFloat)gm_popcount64(oddmask0 | (oddmask1 << 1));
        }
        /*--Adjust for end pad---*/
        const int nfl = vectors->num_field_local;
        const _Bool final_proc = GMEnv_proc_num_field(env) ==
                                 GMEnv_num_proc_field(env)-1;
        const int num_field_inactive_local = final_proc ?
          nfl - (vectors->num_field - vectors->num_field_active) : nfl;
        count -= 2 * num_field_inactive_local;
        /*---Finish---*/
        GMAssert(sum >= 0 && sum <= 2 * vectors->num_field);
        GMAssert(count >= 0 && count <= vectors->num_field);
        sums_local[i] = sum;
        counts_local[i] = count;
      } else { /*---sparse---*/
        for (f = 0; f < vectors->num_packedval_field_local; ++f) {
          /*---Fast way: sum all 64 bits of each word immediately---*/
          const GMBits2x64 value = GMVectors_bits2x64_get(vectors, f, i, env);
          sum += (GMFloat)gm_popcount64(value.data[0]);
          sum += (GMFloat)gm_popcount64(value.data[1]);
        }
        GMAssert(sum >= 0 && sum <= 2 * vectors->num_field);
        sums_local[i] = sum;
      } /*---sparse---*/
    } /*---for i---*/
    /*----------*/
  } /*---if---*/
  /*----------*/

  /*---Do reduction across field procs if needed---*/

  if (num_proc > 1) {
    int mpi_code = 0;
    mpi_code *= 1; /*---Avoid unused variable warning---*/
    mpi_code = MPI_Allreduce(sums_local, sums, vectors->num_vector_local,
                      GM_MPI_FLOAT, MPI_SUM, GMEnv_mpi_comm_field(env));
    GMAssertAlways(mpi_code == MPI_SUCCESS);
    if (env->sparse) {
      mpi_code = MPI_Allreduce(counts_local, counts, vectors->num_vector_local,
                        GM_MPI_FLOAT, MPI_SUM, GMEnv_mpi_comm_field(env));
      GMAssertAlways(mpi_code == MPI_SUCCESS);
    }
  }
}

/*---------------------------------------------------------------------------*/

void GMVectorSums_compute(GMVectorSums* vector_sums,
                          GMVectors* vectors,
                          GMEnv* env) {
  GMAssertAlways(vector_sums != NULL);
  GMAssertAlways(vectors != NULL);
  GMAssertAlways(env != NULL);

  switch (GMEnv_metric_type(env)) {
    /*----------------------------------------*/
    case GM_METRIC_TYPE_SORENSON: {
      /*----------------------------------------*/

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    /*----------------------------------------*/
    case GM_METRIC_TYPE_CZEKANOWSKI: {
      /*----------------------------------------*/
      gm_compute_float_vector_sums(vectors, (GMFloat*)vector_sums->sums,
                                   (GMFloat*)vector_sums->sums_tmp, env);
    } break;
    /*----------------------------------------*/
    case GM_METRIC_TYPE_CCC: {
      /*----------------------------------------*/
      gm_compute_bits2_vector_sums(vectors, (GMFloat*)vector_sums->sums,
                                   (GMFloat*)vector_sums->sums_tmp,
                                   (GMFloat*)vector_sums->counts,
                                   (GMFloat*)vector_sums->counts_tmp, env);
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
