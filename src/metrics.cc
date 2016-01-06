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
#include <math.h>

#include "mpi.h"

#include "env.h"
#include "metrics.h"

#ifdef __cplusplus
extern "C" {
#endif

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

  *metrics = GMMetrics_null();

  if (!Env_is_proc_active(env)) {
    return;
  }

  metrics->data_type_id = data_type_id;
  metrics->num_vector_local = num_vector_local;
  metrics->num_elts_0 = 0;
  metrics->num_elts_01 = 0;

  /*---Compute global values---*/

  const int num_proc = Env_num_proc_vector(env);

  size_t num_vector_bound = 0;
  num_vector_bound =
      num_vector_bound * 1; /*---Avoid unused variable warning---*/
  num_vector_bound = num_proc * (size_t)metrics->num_vector;
  GMAssert(num_vector_bound == (size_t)(int)num_vector_bound
               ? "Vector count too large to store in 32-bit int."
               : 0);

  int mpi_code = 0;
  mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/
  mpi_code = MPI_Allreduce(&(metrics->num_vector_local), &(metrics->num_vector),
                           1, MPI_INT, MPI_SUM, Env_mpi_comm_vector(env));
  GMAssert(mpi_code == MPI_SUCCESS);

  /*---Assume the following to simplify calculations---*/

  GMInsist(
      env,
      num_vector_local >= Env_num_way(env)
          ? "Currently require number of vecs on a proc to be at least num-way"
          : 0);

  const int i_proc = Env_proc_num_vector(env);

  /*---Compute number of elements etc.---*/

  /*--------------------*/
  if (Env_all2all(env)) {
    /*--------------------*/
    /*--------------------*/
    if (Env_num_way(env) == 2) {
      /*--------------------*/
      /*---Store strict upper triang of diag block and half
           the off-diag blocks - a wrapped-rectangle block-row---*/
      metrics->num_elts_local = 0;
      /*---Compute size part 1: (triangle) i_proc==j_proc part---*/
      const int nchoosek = gm_nchoosek(num_vector_local, Env_num_way(env));
      metrics->num_elts_local += nchoosek;
      metrics->num_elts_0 = metrics->num_elts_local;
      /*---Compute size part 2: (wrapped rectangle) i_proc!=j_proc part---*/
      const int num_offdiag_block =
          num_proc % 2 == 0 && 2 * i_proc >= num_proc
              ? (num_proc / 2) - 1
              : (num_proc / 2);
      metrics->num_elts_local +=
          num_offdiag_block * num_vector_local * num_vector_local;
      /*---Allocate index---*/
      metrics->coords_global_from_index =
          (size_t*)malloc(metrics->num_elts_local * sizeof(size_t));
      GMAssert(metrics->coords_global_from_index != NULL);
      /*---Set index part 1: (triangle) i_proc==j_proc part---*/
      size_t index = 0;
      int j = 0;
      for (j = 0; j < num_vector_local; ++j) {
        const size_t j_global = j + num_vector_local * i_proc;
        int i = 0;
        for (i = 0; i < j; ++i) {
          const size_t i_global = i + num_vector_local * i_proc;
          metrics->coords_global_from_index[index++] =
              i_global + metrics->num_vector * j_global;
        }
      }
      /*---Set index part 2: (wrapped rectangle) i_proc!=j_proc part---*/
      const size_t beg = (i_proc + 1) * num_vector_local;
      const size_t end =
          (i_proc + num_offdiag_block + 1) * num_vector_local;
      size_t j_global_unwrapped = 0;
      for (j_global_unwrapped = beg; j_global_unwrapped < end;
           ++j_global_unwrapped) {
        const size_t j_global = j_global_unwrapped % metrics->num_vector;
        int i = 0;
        for (i = 0; i < num_vector_local; ++i) {
          const size_t i_global = i + num_vector_local * i_proc;
          metrics->coords_global_from_index[index++] =
              i_global + metrics->num_vector * j_global;
        }
      }
      GMAssert(index == metrics->num_elts_local);
      /*--------------------*/
    } else /* (Env_num_way(env) == 3) */ {
      /*--------------------*/
      GMAssert(Env_num_way(env) == 3);
      /*---Make the following assumption to greatly simplify calculations---*/
      GMInsist(env, num_proc <= 2 || metrics->num_vector_local % 6 == 0
                        ? "3way all2all case requires num vectors per proc "
                          "divisible by 6."
                        : 0);
      metrics->num_elts_local = 0;
      /*---Compute size pt 1: (tetrahedron) i_proc==j_proc==k_proc part---*/
      const int nchoosek = gm_nchoosek(num_vector_local, Env_num_way(env));
      metrics->num_elts_local += nchoosek;
      metrics->num_elts_0 = metrics->num_elts_local;
      /*---Compute size pt 2: (triang prisms) i_proc!=j_proc==k_proc part---*/
      const int nchoosekm1 = gm_nchoosek(num_vector_local, Env_num_way(env) - 1);
      const int num_procm1 = num_proc - 1;
      metrics->num_elts_local += num_procm1 * nchoosekm1 * num_vector_local;
      metrics->num_elts_01 = metrics->num_elts_local;
      /*---Compute size pt 3: (block sections) i_proc!=j_proc!=k_proc part---*/
      const int num_procm2 = num_proc - 2;
      metrics->num_elts_local += num_procm1 * num_procm2 * num_vector_local *
                                 num_vector_local * (num_vector_local / 6);
      /*---Allocate index---*/
      metrics->coords_global_from_index =
          (size_t*)malloc(metrics->num_elts_local * sizeof(size_t));
      GMAssert(metrics->coords_global_from_index != NULL);
      /*---Set index part 1: (tetrahedron) i_proc==j_proc==k_proc part---*/
      size_t index = 0;
      int k = 0;
      for (k = 0; k < num_vector_local; ++k) {
        const int k_proc = i_proc;
        const size_t k_global = k + num_vector_local * k_proc;
        int j = 0;
        for (j = 0; j < k; ++j) {
          const int j_proc = i_proc;
          const size_t j_global = j + num_vector_local * j_proc;
          int i = 0;
          for (i = 0; i < j; ++i) {
            const size_t i_global = i + num_vector_local * i_proc;
            GMAssert(index < metrics->num_elts_local);
            metrics->coords_global_from_index[index++] =
                i_global +
                metrics->num_vector *
                    (j_global + metrics->num_vector * (k_global));
          }
        }
      }
      /*---Set index part 2: (triang prisms) i_proc!=j_proc==k_proc part---*/
      int k_proc = 0;
      for (k_proc = 0; k_proc < num_proc; ++k_proc) {
        if (k_proc == i_proc) {
          continue;
        }
        const int j_proc = k_proc;
        int k = 0;
        for (k = 0; k < num_vector_local; ++k) {
          const size_t k_global = k + num_vector_local * k_proc;
          int j = 0;
          for (j = 0; j < k; ++j) {
            const size_t j_global = j + num_vector_local * j_proc;
            int i = 0;
            for (i = 0; i < num_vector_local; ++i) {
              const size_t i_global = i + num_vector_local * i_proc;
              GMAssert(index < metrics->num_elts_local);
              metrics->coords_global_from_index[index++] =
                  i_global +
                  metrics->num_vector *
                      (j_global + metrics->num_vector * (k_global));
            }
          }
        }
      } /*---k_proc---*/
      /*---Set index part 3: (block sections) i_proc!=j_proc!=k_proc part---*/
      for (k_proc = 0; k_proc < num_proc; ++k_proc) {
        if (k_proc == i_proc) {
          continue;
        }
        int j_proc = 0;
        for (j_proc = 0; j_proc < num_proc; ++j_proc) {
          if (j_proc == i_proc || j_proc == k_proc) {
            continue;
          }
          const int section_axis = gm_metrics_3way_section_axis(
              metrics, i_proc, j_proc, k_proc, env);
          const int section_num =
              gm_metrics_3way_section_num(metrics, i_proc, j_proc, k_proc, env);
          const int inc = num_vector_local / 6;
          const int k_min = section_axis == 2 ? (section_num)*inc : 0;
          const int k_max =
              section_axis == 2 ? (section_num + 1)*inc : num_vector_local;
          int k = 0;
          for (k = k_min; k < k_max; ++k) {
            const size_t k_global = k + num_vector_local * k_proc;
            const int j_min = section_axis == 1 ? (section_num)*inc : 0;
            const int j_max =
                section_axis == 1 ? (section_num + 1)*inc : num_vector_local;
            int j = 0;
            for (j = j_min; j < j_max; ++j) {
              const size_t j_global = j + num_vector_local * j_proc;
              const int i_min = section_axis == 0 ? (section_num)*inc : 0;
              const int i_max = section_axis == 0 ? (section_num + 1)*inc
                                                  : num_vector_local;
              int i = 0;
              for (i = i_min; i < i_max; ++i) {
                const size_t i_global = i + num_vector_local * i_proc;
                GMAssert(index < metrics->num_elts_local);
                metrics->coords_global_from_index[index++] =
                    i_global +
                    metrics->num_vector *
                        (j_global + metrics->num_vector * (k_global));
              }
            }
          }
        } /*---j_proc---*/
      }   /*---k_proc---*/
      GMAssert(index == metrics->num_elts_local);
    }
    /*--------------------*/
  } else { /*---if not all2all---*/
           /*--------------------*/
    const int nchoosek = num_vector_local >= Env_num_way(env)
                             ? gm_nchoosek(num_vector_local, Env_num_way(env))
                             : 0;
    metrics->num_elts_local = nchoosek;
    metrics->coords_global_from_index =
        (size_t*)malloc(metrics->num_elts_local * sizeof(size_t));
    GMAssert(metrics->coords_global_from_index != NULL);
    /*---LATER: generalize this to N-way---*/
    if (Env_num_way(env) == 2) {
      /*---Need store only strict upper triangular part of matrix---*/
      size_t index = 0;
      int j = 0;
      for (j = 0; j < num_vector_local; ++j) {
        const size_t j_global = j + num_vector_local * i_proc;
        int i = 0;
        for (i = 0; i < j; ++i) {
          const size_t i_global = i + num_vector_local * i_proc;
          metrics->coords_global_from_index[index++] =
              i_global + metrics->num_vector * j_global;
        }
      }
      GMAssert(index == metrics->num_elts_local);
    } else /* (Env_num_way(env) == 3) */ {
      GMAssert(Env_num_way(env) == 3);
      /*---Need store only strict interior of tetrahedron---*/
      size_t index = 0;
      int k = 0;
      for (k = 0; k < num_vector_local; ++k) {
        const size_t k_global = k + num_vector_local * i_proc;
        int j = 0;
        for (j = 0; j < k; ++j) {
          const size_t j_global = j + num_vector_local * i_proc;
          int i = 0;
          for (i = 0; i < j; ++i) {
            const size_t i_global = i + num_vector_local * i_proc;
            GMAssert(index < metrics->num_elts_local);
            metrics->coords_global_from_index[index++] =
                i_global +
                metrics->num_vector *
                    (j_global + metrics->num_vector * (k_global));
          }
        }
      }
      GMAssert(index == metrics->num_elts_local);
    } /*---if num_way---*/
    /*--------------------*/
  } /*---if all2all---*/
  /*--------------------*/

  /*---Allocations---*/

  switch (data_type_id) {
    /*----------*/
    case GM_DATA_TYPE_BITS1: {
      /*---(design not complete)---*/
      const size_t num_floats_needed =
          gm_ceil_i8(metrics->num_elts_local, 8 * sizeof(GMFloat));
      metrics->data = malloc(num_floats_needed * sizeof(GMFloat));
      GMAssert(metrics->data != NULL);
      metrics->data_type_num_values = 1;
    } break;
    /*----------*/
    case GM_DATA_TYPE_FLOAT:
//---TODO: remove thee casts on the mallocs, here and elsewhere.
      metrics->data = (GMFloat*)malloc(metrics->num_elts_local*sizeof(GMFloat));
      GMAssert(metrics->data != NULL);
      metrics->data_type_num_values = 1;
      break;
    /*----------*/
    case GM_DATA_TYPE_TALLY2X2: {
      metrics->data
        = (GMTally2x2*)malloc(metrics->num_elts_local*sizeof(GMTally2x2));
      GMAssert(metrics->data != NULL);
      metrics->data_M
        = (GMFloat*)malloc(metrics->num_elts_local*sizeof(GMFloat));
      GMAssert(metrics->data_M != NULL);
      metrics->data_type_num_values = 4;
    } break;
    /*----------*/
    case GM_DATA_TYPE_TALLY4X2: {
      metrics->data
        = (GMTally4x2*)malloc(metrics->num_elts_local*sizeof(GMTally4x2));
      GMAssert(metrics->data != NULL);
      metrics->data_M
        = (GMFloat*)malloc(metrics->num_elts_local*sizeof(GMFloat));
      GMAssert(metrics->data_M != NULL);
      metrics->data_type_num_values = 8;
    } break;
    /*----------*/
    default:
      GMAssert(GM_BOOL_FALSE ? "Invalid data type." : 0);
  } /*---switch---*/
}

/*===========================================================================*/
/*---Metrics pseudo-destructor---*/

void GMMetrics_destroy(GMMetrics* metrics, GMEnv* env) {
  GMAssert(metrics);
  GMAssert(env);
  GMAssert(metrics->data != NULL || !Env_is_proc_active(env));

  if (!Env_is_proc_active(env)) {
    return;
  }

  free(metrics->data);
  free(metrics->coords_global_from_index);
  if (metrics->data_M != NULL) {
    free(metrics->data_M);
  }
  *metrics = GMMetrics_null();
}

/*===========================================================================*/
/*---Metrics checksum---*/

/*---Helper function - perform one bubble sort step---*/

static void gm_bubbledown(size_t* i, size_t* j) {
  if (*i < *j) {
    const size_t tmp = *i;
    *i = *j;
    *j = tmp;
  }
}

/*---------------------------------------------------------------------------*/

static size_t gm_lshift(size_t a, int j) {
  if (j >= 64  || j <= -64) {
    return 0;
  }
  return j > 0 ? a << j : a >> (-j);
}

/*---------------------------------------------------------------------------*/

GMChecksum GMMetrics_checksum(GMMetrics* metrics, GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(metrics->data != NULL || !Env_is_proc_active(env));

  /*---Initializations---*/

  GMChecksum result = GMChecksum_null();

  if (!Env_is_proc_active(env)) {
    return result;
  }

  GMAssert(Env_num_way(env) <= 3 ? "This num_way not supported." : 0);

  /*---Initializations---*/

  typedef size_t UI64;
  GMStaticAssert(sizeof(UI64) == 8);

  UI64 sums_l[16];
  int i = 0;
  for (i = 0; i < 16; ++i) {
    sums_l[i] = 0;
  }
  double sum_d = 0;

  const int w = 30;
  GMAssertAlways(64-2*w >= 4);
  const UI64 one64 = 1;

  const UI64 lomask = ( one64 << w ) - 1;
  const UI64 lohimask = ( one64 << (2*w) ) - 1;

  UI64 coords[3];
  UI64 index = 0;
  for (index = 0; index < metrics->num_elts_local; ++index) {
    /*---Obtain global coords of metrics elt---*/
    coords[2] = 0;
    for (i = 0; i < Env_num_way(env); ++i) {
      coords[i] =
          GMMetrics_coord_global_from_index(metrics, index, i, env);
    }
    /*---Reflect coords by symmetry to get uniform result---*/
    gm_bubbledown(&coords[1], &coords[2]);
    gm_bubbledown(&coords[0], &coords[1]);
    gm_bubbledown(&coords[1], &coords[2]);

    /*---Loop opver data values at this index---*/
    int i_value = 0;
    for (i_value=0; i_value<metrics->data_type_num_values; ++i_value) {
        /*---Construct global id for metrics data vbalue---*/
        UI64 uid = coords[0];
        for (i = 1; i < Env_num_way(env); ++i) {
          uid = uid * metrics->num_vector + coords[i];
        }
        uid = uid * metrics->data_type_num_values + i_value
        /*---Randomize---*/
        const UI64 rand1 = gm_randomize(uid + 956158765);
        const UI64 rand2 = gm_randomize(uid + 842467637);
        UI64 rand_value = rand1 + gm_randomize_max() * rand2;
        rand_value &= lohimask;
        /*---Now pick up value of this metrics elt---*/
        GMFloat value = 0;
        switch (metrics->data_type_id) {
          /*--------------------*/
          case GM_DATA_TYPE_BITS1: {
            GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
          } break;
          /*--------------------*/
          case GM_DATA_TYPE_FLOAT: {
            value = GMMetrics_czekanowski_get_from_index(metrics, index, env);
          } break;
          /*--------------------*/
          case GM_DATA_TYPE_TALLY2X2: {
            const int i0 =  i_value / 2;
            const int i1 =  i_value % 2;
            value = GMMetrics_ccc_get_from_index_2(metrics, index, i0, i1, env);
          } break;
          /*--------------------*/
          case GM_DATA_TYPE_TALLY4X2: {
            const int i0 =  i_value / 4;
            const int i1 = (i_value / 2 ) % 2;
            const int i2 =  i_value % 2;
            value = GMMetrics_ccc_get_from_index_3(metrics, index,
                                                   i0, i1, i2, env);
          } break;
          /*--------------------*/
          default:
            GMAssertAlways(GM_BOOL_FALSE ? "Invalid data type." : 0);
        } /*---switch---*/
        /*---Convert to integer.  Store only 2*w bits max---*/
        const int log2_value_max = 4;
        GMAssertAlways(value >= 0 && value < (1<<log2_value_max));
        UI64 ivalue = value * (one64 << (2*w-log2_value_max));
        /*---Multiply the two values---*/
        const UI64 a = rand_value;
        const UI64 alo = a & lomask;
        const UI64 ahi = a >> w;
        const UI64 b = ivalue;
        const UI64 blo = b & lomask;
        const UI64 bhi = b >> w;
        const UI64 cx = alo * bhi + ahi * blo;
        UI64 clo = alo * blo + ( (cx & lomask) << w);
        UI64 chi = ahi * bhi + (cx >> w);
        sum_d += ivalue * (double)rand_value / ((double)(one64<<(2*w)));
        /*---(move the carry bits)---*/
        chi += clo >> (2*w);
        clo &= lohimask;
        /*---Split the product into one-char chunks, accumulate to sums---*/
        for (i = 0; i < 8; ++i) {
          sums_l[0+i] += (clo << (64-8-8*i)) >> (64-8);
          sums_l[8+i] += (chi << (64-8-8*i)) >> (64-8);
        }
    } /*---for i_value---*/
  } /*---for index---*/
  /*---Global sum---*/
  UI64 sums_g[16];
  int mpi_code = 0;
  mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/
  mpi_code = MPI_Allreduce(sums_l, sums_g, 16, MPI_UNSIGNED_LONG_LONG,
                           MPI_SUM, Env_mpi_comm_vector(env));
  GMAssertAlways(mpi_code == MPI_SUCCESS);

  /*---Combine results---*/

  for (i = 0; i < GM_CHECKSUM_SIZE; ++i) {
    int j = 0;
    for (j = 0; j < 8; ++j) {
      result.data[i] += gm_lshift(sums_g[0+j], 8*j-2*w*i) & lohimask;
      result.data[i] += gm_lshift(sums_g[8+j], 8*j-2*w*(i-1)) & lohimask;
    }
  }
  /*---(move the carry bits---*/
  result.data[1] += result.data[0] >> (2*w);
  result.data[0] &= lohimask;
  result.data[2] += result.data[1] >> (2*w);
  result.data[1] &= lohimask;

  /*---Check against floating point result---*/

  const double tmp = sum_d;
  mpi_code = MPI_Allreduce(&tmp, &sum_d, 1, MPI_DOUBLE, MPI_SUM,
                           Env_mpi_comm_vector(env));
  GMAssertAlways(mpi_code == MPI_SUCCESS);

  double result_d = result.data[0] / ((double)(one64<<(2*w))) +
                    result.data[1] +
                    result.data[2] * ((double)(one64<<(2*w)));
  result_d = 1 * result_d; /*---Avoid unused variable warning---*/
  GMAssertAlways(fabs(sum_d-result_d) < sum_d * 1.e-10);

  return result;
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
